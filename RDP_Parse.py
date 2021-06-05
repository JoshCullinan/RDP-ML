import numpy
import re
import pandas as pd
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Parse RDP Output from SantaSim')
    parser.add_argument('--f', dest='filename', type=str,help='File to Parse')
    args = parser.parse_args()
    
    #Create relative path names
    import os
    dirname = os.path.dirname(__file__)
    
    filename = 'alignment_3.fa.csv'

    RDP_Output = pd.read_csv(os.path.join(dirname, filename), delimiter='\n')

    for index, row in Events.iterrows():
        out = row.item()[row.item().find("[")+1:row.item().find("]")]
        if out != '':
            Events.iloc(axis=0)[index] = out
        else:
            Events.iloc(axis=0)[index] = None

    Events = Events.squeeze()

    Events = Events.str.split(pat = ', ', expand = True)

    Events.to_csv(path_or_buf = filename + '_Parsed_.csv', sep = ',')

