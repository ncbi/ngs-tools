************** STANDARD OUT FORMAT *******************
The default out format for the aligns with a dbs version looks like number of lines like this:
1026    10239
1030    9526x2
1031    10376x3
1115    10376   314293x2
1149    10239   10376

Where all values are tab-separated. The first number is a spot id, then goes list of tax ids.
Here:
spot 1026 has one kmer hit to 10239
spot 1030 has 2 kmers that belongs to tax id 9526
spot 1115 has 1 kmer hit to 10376 and 2 kmers of 314293
and so on

Sometimes number of hits are not important so aligns_to support optional [-hide_counts] flag and with this flag provided the output will look like this:
1026    10239
1030    9526
1031    10376
1115    10376   314293
1149    10239   10376

Basically, it's same, but without count information.
Spot numbers extracted from the SRA object directly, so they are integer numbers, but when working with fasta parsing spot id from the 
sequence description is not a trivial problem, so for fasta files spot id is everything before space or '/', so output file may look like this:

SRR9694553.1026    10239
SRR9694553.1030    9526
SRR9694553.1031    10376
SRR9694553.1115    10376   314293
SRR9694553.1149    10239   10376

WARNING: In fact, some SRA runs use GUIDs as spot id when dumped as fasta objects, so it's wrong to make an assumption that spot id is any integer number.
WARNING: Due to the multithreaded nature of the aligns_to executable it can report spots in any order, it's not going to match the order of the SRA object or any input fasta file.

**************** COMPACT FORMAT ****************
Large runs can produce very long out files, so the compact out format was introduced. It's supported with [-compact] optional flag for the aligns_to executable.
It's a little bit more complicated than the default format and looks like this:
(look for "void print_compact" function for the source code)
1       2759
1       9443
5       9526
1       9604
3       9606
18      10239
1       10239   10376
1       10239   2364194
6       10376
1       11632
1       11676
2       131567
2       207598
2       314293
1       314295
1       2759
1       9347
1       9443
3       9526
1       9606
15      10239
1       10239   2364194
4       10376
1       32524
1       131567
4       207598
2       314146
4       314293
2       314295
        5632    10239
        5632    2365021
3       2759
1       9443
4       9526
19      10239
4       10376

Input data processed in blocks and information for each block reported separately one next to another. 
15      10239 means that this block has 15 spots with tax id 10239
1       10239   2364194  means that there is 1 spot with this set of tax ids. This spot is not included into the above spot summary

So with the above 2 lines there are 16 spots total identified and one of them has hits to both (10239   2364194) and 15 others have hits to 10239 only.

WARNING: Don't make any assumptions on the block size, this can change in future.
The main complexity of this format is that it also may have such lines as
        5632    10239
        5632    2365021

They start with the tab character. This means that aligns_to failed to merge such reads into spots (because they may belong to different blocks) and 
such lines reported in the different format (spot id, set of tax ids) described at the beginning of the readme file.

