import pandas as pd 
import numpy as np
import ast
import argparse
from collections import defaultdict
from pathlib import Path
import os

class rdpParser():
	rdpcsv_fileName = ""
	seqmap_fileName = ""
	events_fileName = ""
	maximum_genome_length = 3838

	seqmap_df = pd.DataFrame()
	events_df = pd.DataFrame()
	rdp_df = pd.DataFrame()
	samplesize = 0
	added_events_count = 0
	seqmap_dict = {}
	events_dict = {}
	events_map = {}
	inv_seqmap_dict = defaultdict(set)
	

	#variables to test rdp accuracy
	true_recombinants = []
	total_recombinants = 0
	matched_recombinants = 0	
	

	def __init__(self, rdpf, seqmapf, eventsf):
		self.rdpcsv_fileName = rdpf
		self.seqmap_fileName = seqmapf
		self.events_fileName = eventsf
		self.read_files()
		self.create_dictionaries()
		self.add_extra_events()
	
	#reads csv files into pandas dataframes
	def read_files(self):
		self.rdp_df = pd.read_csv(self.rdpcsv_fileName, index_col=False)
		#remove whitespaces because rdp csv has spaces before names sometimes
		self.rdp_df = self.rdp_df.rename(columns=lambda x: x.strip())
		self.seqmap_df = pd.read_csv(self.seqmap_fileName, delimiter='*', index_col='Sequence')
		self.events_df = pd.read_csv(self.events_fileName, delimiter='*', usecols= ['EventNum','Breakpoints'])
		self.samplesize = self.seqmap_df.tail(1).index.item()

	#constructs dictionaries from dataframes
	def create_dictionaries(self):
		#generating dictionaries from dataframes 
		self.seqmap_dict = {i:ast.literal_eval(v) for i, v in enumerate(self.seqmap_df['Events'].to_numpy(), 1)}
		self.events_dict = {event:ast.literal_eval(bp) for event,bp in zip(self.events_df['EventNum'], self.events_df['Breakpoints'])}

		#Creating an inverted seqmap dictionary (event:sequences instead of sequence:events) 
		#with key,value pairs: event, [sequences containing event]
		for key, value in self.seqmap_dict.items():
			for eventnum in value:
				self.inv_seqmap_dict[eventnum].add(key)

	#adds extra events to account for events present in majority of genomes
	#also constructs the final dictionary which is used to compare sim data with rdp output
	def add_extra_events(self):
		#Checking this dictionary for events in >50% of sequences. For these we will create a new event which is identical, 
		#except present in the complement set of sequences.

		if (len(self.inv_seqmap_dict.keys()) > 0):
			max_event_num = max(self.inv_seqmap_dict.keys())
		else:
			max_event_num = 0
		
		counter = 0
		for key, value in self.inv_seqmap_dict.items():
			if (len(value) >= self.samplesize/2):
				counter += 1
				#adding new event to events dictionary
				self.events_dict[max_event_num+counter] = self.events_dict[key]

				#adding new event to event:sequence dictionary
				complemented_sequences = set(range(1,self.samplesize+1)) - value	
				for seq in complemented_sequences:
					self.inv_seqmap_dict[max_event_num+counter].add(seq)

		self.added_events_count = counter;

		for k,v in self.events_dict.items():
			self.events_map[k] = [v, self.inv_seqmap_dict[k]]


	def compare_rdp_with_sim(self):
		
		if (len(self.rdp_df) > 0):
			num_of_events = self.rdp_df.at[self.rdp_df.index[-1], 'Event']
		else:
			num_of_events = 0	
	
		for i in range(num_of_events):
			#collect breakpoints from rdp output
			#collect sequence sets for i,i+1,i+2 from rdp output
			#check events map for events that match bp, and contain recombinants inside one of the 3 sequence sets
			#the sequence set that contains those recombinants is the recombinant sequence

			#if there is a match found, add it to output file. Also note how many of the recombinants were identified			
			#if there is no match inside events map, list it as false positive and dont add it to outputfile
			#also keep track of how many events were correctly identified, taking into account added_events_count

			#collecting breakpoints and sequence sets from RDP output
			y = i*3
			start_bp = self.rdp_df.at[y, 'StartBP'] 
			end_bp = self.rdp_df.at[y, 'EndBP']	
			#3 sets of sequences for recombinant+parents triplet
			seqs = [{int(j) for j in self.rdp_df.at[y+k, 'ISeqs(A)'].split('$')} for k in range(3)]

			#swapping bp if "wrong" way around (circular genome)
			if end_bp < start_bp:
				start_bp, end_bp = end_bp, start_bp
			
			#now search simulated events for matching breakpoints and sequences	
			match_found = False		
			for eventnum, bp_recombinant_list in self.events_map.items():
				breakpoints = bp_recombinant_list[0]
				recombinants = bp_recombinant_list[1]

				#also need to consider for some events RDP will wrap around the end (circular genome)
				#while being very near the end of the genome
				#thus using endpoint instead of 0
				#second condition accounts for this happening.
				first_condition = (breakpoints[0]-100 <= start_bp and start_bp <= breakpoints[0]+100) and (breakpoints[1]-100 <= end_bp and end_bp <= breakpoints[1]+100)
				second_condition = False				

				if (self.maximum_genome_length-100 <= end_bp and end_bp <= self.maximum_genome_length+100):
					second_condition = (breakpoints[0]-100 <= 0 and 0 <= breakpoints[0]+100) and (breakpoints[1]-100 <= start_bp and start_bp <= breakpoints[1]+100)
			

				if (first_condition or second_condition):
					for index, seq_set in enumerate(seqs,0):
						if not recombinants.isdisjoint(seq_set):
							#match found
							self.total_recombinants += len(recombinants)
							self.matched_recombinants += len(recombinants.intersection(seq_set))

							#true_recombinants lists all matched events (y = index of event in RDP output),
							#as well as which of the 3 candidate sequence sets are the true recombinant
							self.true_recombinants.append([y, index])							
							break

	def export_ml_data(self):			

		self.rdp_df = self.rdp_df.drop(['Event', 'StartBP', 'EndBP', 'ISeqs(A)'], axis=1)
		file_name = "output/ml_input.txt"
		file_path = Path(file_name)

		#if file doesnt exist yet, create it and write header
		if (not (file_path.exists())):	
			os.makedirs(os.path.dirname(file_name), exist_ok=True)		
			with open(file_name, "w+") as g:
				header = ['Recombinant'] + self.rdp_df.columns.to_list()
				g.write('\t'.join(str(s) for s in header) + '\n')


		with open(file_name, "a") as f:
			for y, r in self.true_recombinants:
				out = [self.rdp_df.iloc[y+s].tolist() for s in range(3)]	
				zipped = zip(out[0], out[1], out[2])

				recombinant_vector = [0.0,0.0,0.0]
				recombinant_vector[r] = 1.0
				final_output = [tuple(recombinant_vector)] + [x for x in zipped]

				f.write('\t'.join(str(y) for y in final_output) + '\n')

	def export_rdp_accuracy_data(self, rdp_recIDTests_filename):
		out = []

		total_event_number = len(self.events_map.keys()) - self.added_events_count
		matched_events = len(self.true_recombinants)

		out.append(total_event_number)
		out.append(matched_events)		

		#of the matched recombinant events, in how many of the triplets was the recombinant correctly identified?
		#retrieving consensus values from rdp output, and then for each matched recombinant, which one of the triplet is max?		
		rdp_rec_df = pd.read_csv(rdp_recIDTests_filename, index_col=False)		
		rdp_rec_df = rdp_rec_df.rename(columns=lambda x: x.strip())

		#find which triplet has highest consensus (0,1 or 2), for all events that have been matched
		rdp_consensus = [ rdp_rec_df.iloc[y:y+3, 3].idxmax(axis=0)-y for y in [x for x,y in self.true_recombinants]]

		#now compare these to the simulation results
		correctly_identified_events = 0
		for i, r in enumerate(rdp_consensus):
			if (r == self.true_recombinants[i][1]):
				correctly_identified_events += 1

		out.append(correctly_identified_events)

		#total amount of recombinants from events identified by RDP (including multiple from one event)
		#and how many of these were identified by RDP
		out.append(self.total_recombinants)
		out.append(self.matched_recombinants)

		#check if file exists, if not, create file + header
		file_name = "output/rdp_accuracy.txt"
		file_path = Path(file_name)

		#if file doesnt exist yet, create it and write header
		if (not (file_path.exists())):
			os.makedirs(os.path.dirname(file_name), exist_ok=True)			
			with open(file_name, "w+") as g:
				header = ['TotalEvents', 'MatchedEvents', 'RecombinantsCorrectlyChosen', 'TotalRecombinants', 'MatchedRecombinants']
				g.write('\t'.join(str(s) for s in header) + '\n')

		with open(file_name, "a") as file:
			file.write('\t'.join(str(y) for y in out) + '\n')		


#Parsers command line input to get file names. Pipeline.bat currently inputs filenames. Can be done manually as well.
def getFileNames():
	#Define parser
	parser = argparse.ArgumentParser(description='Parse Recombination Information from SantaSim')
	
	#Add arguments for command line
	parser.add_argument('--rdpcsv', dest = 'rdpcsv_fileName', type=str, help='RDP Recombination Identifying Stats File', required=True)
	parser.add_argument('--seq', dest = 'seqmap_fileName', type=str,  help='Sequence Events file from Simulation', required=True)
	parser.add_argument('--rec', dest = 'events_fileName', type=str,  help='Recombinaton events file from Simulation', required=True)
	parser.add_argument('--IDTests', dest = 'rdp_recIDTests_fileName', type=str,  help='RDP ID Tests file', required=True)
	
	#Parse Events
	args = parser.parse_args()

	#Returns parsed file names. Could be done in a 'nicer' manner but it works for now.
	return(args.rdpcsv_fileName, args.seqmap_fileName, args.events_fileName, args.rdp_recIDTests_fileName)
			
if __name__=='__main__':
	rdpcsv_fileName,seqmap_fileName, events_fileName, rdp_recIDTests_fileName = getFileNames()

	parser = rdpParser(rdpcsv_fileName, seqmap_fileName, events_fileName)
	parser.compare_rdp_with_sim()
	parser.export_ml_data()
	parser.export_rdp_accuracy_data(rdp_recIDTests_fileName)


### If needed for testing or if parser fails ###
	# rdpcsv_fileName = "alignment_1.faRecombIdentifyStats.csv"
	# seqmap_fileName = "sequence_events_map_1.txt"
	# events_fileName = "recombination_events_1.txt"
	# rdp_recIDTests_fileName = "alignment_1.faRecIDTests.csv"