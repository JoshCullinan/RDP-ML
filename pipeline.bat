@REM Set EPOCH to the number of outputs required
set /A EPOCH=3

call conda activate RDP

@REM muscle3.8.31_i86win32.exe -in alignment_1.fa -out MuscleAligned.afa
FOR /L %%A IN (2,1,%EPOCH%) DO (

    java -jar santa.jar low_recomb_rate.xml

    Ren "alignment_1.fa" "alignment_%%A.fa"
    Ren "recombination_events.txt" "recombination_events_%%A.txt"
    Ren "sequence_events_map.txt" "sequence_events_map_%%A.txt"
    
    python Parse_Events.py --f sequence_events_map_%%A.txt

    RDP5CL.exe -f alignment_%%A.fa REM -am -o -nor

)
