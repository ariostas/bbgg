Here you can find all my code used for the bbgg analysis.

The main file is the one named bbgg.C, which runs over the samples given as a parameter and using the category also given as parameter. It runs over the small ntuples saved in /afs/cern.ch/work/a/ariostas/public/bbgg/ . The small ntuples were created with the code found on the Selection folder. Available categories are "HighPt" and "LowPt".

To run use root -b -q bbgg.C+ to run over all samples using the high-pt category or use root -b -q bbgg.C+\(\"name of background sample\",\"category\"\) to specify samples and category.

In Results.pdf you can find the event yields and an explanation of what I've done.
