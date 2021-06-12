@echo OFF
@REM Set EPOCH to the number of outputs required
set /A EPOCH=1

@REM Creating a "unique" identifier based on the current time for the alignment files. If you know how to do it better please do - this code is crap. I don't know BAT.
SET T=%date%-%time%
for /f "tokens=1-3 delims=/" %%I in ("%T%") do @set T=%%I-%%J-%%K
for /f "tokens=1-3 delims=:" %%I in ("%T%") do @set T=%%I-%%J-%%K
for /f "tokens=1-2 delims=." %%I in ("%T%") do @set T=%%I-%%J
call conda activate RDP

FOR /L %%A IN (1,1,%EPOCH%) DO (

    @REM Run simulation.
    java -jar santa.jar low_recomb_rate.xml

    @REM Housekeeping to allow for easy file management.
    Ren "alignment_0.fa" "alignment_%%A.fa"
    Ren "recombination_events.txt" "recombination_events_%%A.txt"
    Ren "sequence_events_map.txt" "sequence_events_map_%%A.txt" 
    
    @REM Run RDP scan on the simulation
    RDP\RDP5CL.exe -f ../alignment_%%A.fa -ds

    @REM Parse the output files
    python output_parser.py --rdpcsv "alignment_%%A.faRecombIdentifyStats.csv" --seq "sequence_events_map_%%A.txt" --rec "recombination_events_%%A.txt" --IDTests "alignment_%%A.faRecIDTests.csv" 

    @REM
    md "output\alignments" > nul 
    md "output\alignments\%T%" > nul

    @REM If we decide to delete these file we can just uncomment this line. 
    @REM del "alignment_%%A.faRecombIdentifyStats.csv", "alignment_%%A.faRecIDTests.csv", "alignment_%%A.fa.csv", "alignment_%%A.fa.rdp5"

    @REM Move all of the output files generated into the output folder.
    move "*.csv" "..\RDP-ML\output\alignments\%T%" > nul
    move "*.txt" "..\RDP-ML\output\alignments\%T%" > nul
    move "*.fa" "..\RDP-ML\output\alignments\%T%" > nul
    move "*.rdp5" "..\RDP-ML\output\alignments\%T%" > nul
)

@REM Run the command below to test the output parser from your terminal/CMD, assuming that you've run the above code and have some RDP files to parse.  
    @REM python output_parser.py --rdpcsv "alignment_1.faRecombIdentifyStats.csv" --seq "sequence_events_map_1.txt" --rec "recombination_events_1.txt" --IDTests "alignment_1.faRecIDTests.csv"