@REM Set EPOCH to the number of outputs required
set /A EPOCH=1

call conda activate RDP

FOR /L %%A IN (1,1,%EPOCH%) DO (

    java -jar santa.jar low_recomb_rate.xml

    Ren "alignment_0.fa" "alignment_%%A.fa"
    Ren "recombination_events.txt" "recombination_events_%%A.txt"
    Ren "sequence_events_map.txt" "sequence_events_map_%%A.txt" 
    
    RDP5CL.exe -f alignment_%%A.fa -ds

    python output_parser.py --rdpcsv "alignment_%%A.faRecombIdentifyStats.csv" --seq "sequence_events_map_%%A.txt" --rec "recombination_events_%%A.txt" --IDTests "alignment_%%A.faRecIDTests.csv" 
)

@REM Run the command below to test the output parser from your terminal/CMD, assuming that you've run the above code and have some RDP files to parse.  
    @REM python output_parser.py --rdpcsv "alignment_1.faRecombIdentifyStats.csv" --seq "sequence_events_map_1.txt" --rec "recombination_events_1.txt" --IDTests "alignment_1.faRecIDTests.csv" 
