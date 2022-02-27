# batch_4DTV_calculation
**Calculate 4DTV (transversion rate on 4-fold degenerated sites) much more faster with batch axt files**

This script was based on the "[calculate_4DTV_correction.pl](https://github.com/JinfengChen/Scripts/blob/master/FFgenome/03.evolution/distance_kaks_4dtv/bin/calculate_4DTV_correction.pl)"

The old script can only calculate 4DTV for a pair of sequences at a time, which are contained in an axt file.
In general, there is nothing wrong with this approach.
However, in the case of heavy computation, this not only slows down the progress, but also tends to cause the computer to crash

## Merge axt files
If you have many axt files:

example1.axt
```
seq1-seq2
ATGTCTCATATGTCTTCTGTGAACGCGAAAAATCTTCAAAAACTAGCAGATTCAATTGTC
AAACATGTAAAGCACTTTAACAATAATGAAGTTTTGTGTCTGATCAAACTCTTCAATGTG
CTGATGGGAGAGCAGAGCGAGCACAGGGTTGGAAATGGACTGGATCGTGGTAAATTCAGG
AGCATCCTCCACAACACATTTGGAATGACAGATGACATGATTATGGACAGAGTCTTCCGT
GCATTTGACAAGGACAATGATAGCAACGTCAGTGTAAAAGAATGGATAGAAGGACTTTCA
GTGTTTCTGCGAGGGACCTTGgatgaaaaaattaaatATTGTTTTGAGGTTTATGACTTA
AATGGGGATGGATATATTTCACGAGAAGAGATGTTCCACATGCTGAAAAACAGTCTCATA
AAACAACCAACAGAAGAAGATCCAGATGAAGGGATTAAGGACTTGGTAGAGATTACTCTT
AAAAAGATGgaCCACGATCACGACAGCAGACTTTCATACGCTGATTTTGAGAAAGCAGTA
AAAGAAGAAAATCTCTTGCTTGAGGCTTTTGGAGCTTGTCTTCCTGATGCAAAGagtaTT
CTTGCTTTTGAAAGACAGGCCTTCCAG------GATACCACAGAAAAT
atgctgaaaatgTCGGCGATGAACAGAAAATTAATTCAAAACCTCGCCGAGACTTTATGC
AGACAAGTCAAACATTTTAATAAAACAGAGACGGAGTGTCTGATAAGGCTGTTCAACAGT
CTGCTGGGAGAGCAGGCAGAGAGAAAGACGACTATTGGAGTGGACCGGGCCAAATTCAGA
AATATACTGCACCACACTTTCGGGATGACCGACGACATGATGACGGACAGAGTTTGTCGT
GTCATTGACAAGGACAACGATGGCTACTTAAGCGTTAAAGAGTGGGTTGAggctctgtct
gtctttctaagAGGCACACTGGATGAAAAAATGAAATaCTGTTTTGAGGTGTATGACCTG
AACGGGGATGGATACATCTCACGTGAGGAGATGTTTCAGATGCTGAAAGACAGCCTCATC
AGGCAGCCCACCGAAGAGGATCCTGATGAGGGGATTAAGGATATTGTGGAGATTGCCTTG
AAAAAAATGGATTATGACCATGATGGAAGAGTTTCTTATGCTGATTTTGAGAAGACGGTC
ATGGATGAAAACCTTTTACTAGAAGCTTTTGGAAACTGCCTTCCTGATGCAAAGAGTGTA
CTAGCATTTGAGCAACAGGCATTCCAGAAACACGAACACTGCAAAGAA
```
example2.axt
```
seq3-seq4
ATGGATCGCCATTCCAAtttaatttccatttggctgcaACTGGAACTGTGTGCCATGGCA
GTACTTCTGGCAAAAGGGGAGATAAGATGCTACTGTGATGCAGCGCATTGTGTGGCAACA
GGTTACATGTGTAAATCCGAGTTAAATGCCTGCTTCACCAGGCTTCTGGACCCACAGAAC
ACAAACTCCCCTCTCACGCATGGCTGCTTGGACCCGACTGCAAACACAGCAGATGTTTGC
CATGCTGGAAGGACAGAGAGCCGCGCTGGGGCCTCGGAGAAGCTTGAGTGCTGTCACGAC
GATATGTGCAATTACAGAGGACTCCATGATGTTGTTTCATATCCCAGGGGGGACAGCTCA
GATCATGGAACAAGATATCAGCCAGACAGTAGCAGGAATCTTCTGACCAGGGTTCAGGAT
TTAACATCCTCTAAAGAGCTGTGGTTCAGAGCAGCCGTGATCGCTGTGCCCATCGCTGGG
GGGCTCATTCTAGTGCTTCTCATCATGCTCGCCTTGCGGATGCTTCGAAGTGAAAACAAA
AGACTGCAGGACCAGAGGCAGCAGATGCTGTCCCGCTTGCACTACAACTTTCATGGA---
CACCACACGAAGAAGGGCCAGGTAGCCAAACTGGATTTGGAATGCATGGTTCCCGTAACC
GGACACGAGAACTGCTGTATGACTTGCGACAAACTGCGACAGTCTGAACTCCACAAT---
---------------GATAAATTGCTGTCTTTAGTTCACTGGGGAATTTACAGCGGTCAC
GGGAAATTGGAATttgta
ATGGATCGC---------CTGGTTTCTCTGTGGTTTCAGCTGGAACTTTGTGCGATGGCT
GTTCTTCTCACGAAAGGAGAGATCAGGTGCTACTGTGACGCACCGCACTGCGTTGCCACC
GGATACATGTGTAAATCAGAGCTCAACGCTTGCTTTACTAAGGTCCTGGACCCTCTTAAC
ACAAACTCACCTTTAACACACGGCTGCGTGGATTCGCTTTTAAACTCTGCAGACGTGTGC
TCTAGTAAAAATGTGGACATTTCAAGTGGAAGCTCCTCTCCTGTGGAGTGCTGCCATGAT
GATATGTGTAACTACAGGGGTTTGCATGAC---CTCACACACCCCAGAGGGGACTCAACA
GAC---------CGATACCACAGC---TCCAATCAGAACCTGATCACAAGGGTGCAAGAG
TTAGCGTCTGCTAAAGAGGTGTGGTTCCGGGCGGCGGTGATAGCGGTTCCCATCGCGGGT
GGGCTTATCCTGGTTCTGCTGATTATGCTGGCGTTGCGAATGCTCCGTAGCGAAAACAAG
CGTCTCCAGGCACAGCGCCAGCAGATGCTTTCTCGCCTGCATTACAGCTTTCACGGACAC
CACCATGCCAAGAAAGGCCACGTGGCTAAGTTGGACTTGGAGTGTATGGTGCCGGTAACG
GGACATGAGAACTGTTGTCTGGGCTGCGATAAGCTGCGGCAGACGGATTTGTGCACTGGA
GGAGGAAGCGGGGGTGAGCGTCTCCTATCTCTGGTACACTGGGGGATGTACACGGGGCAC
GGAAAGCTGGAGTTCGTA
```
...

To batch calculate 4DTV, simply merge many axt files into one file (AXT file) using a shell script. 
Note: A sequence does not have a line break in the merged AXT file.

```
> Merged.AXT
for file in `ls  *.axt`;do
  seq=`sed '1d' $file | tr -d "\n"`
  n=`echo $seq | awk '{print length()/2}'`
  m=$((2 * n))
  head -1 $file >> Merged.AXT
  echo ${seq:0:$n} >> Merged.AXT
  echo ${seq:$n:$m} >> Merged.AXT
done
```

Merged.AXT
```
seq1-seq2
ATGTCTCATATGTCTTCTGTGAACGCGAAAAATCTTCAAAAACTAGCAGATTCAATTGTC AAACATGTAAAGCACTTTAACAATAATGAAGTTTTGTGTCTGATCAAACTCTTCAATGTG CTGATGGGAGAGCAGAGCGAGCACAGGGTTGGAAATGGACTGGATCGTGGTAAATTCAGG AGCATCCTCCACAACACATTTGGAATGACAGATGACATGATTATGGACAGAGTCTTCCGT GCATTTGACAAGGACAATGATAGCAACGTCAGTGTAAAAGAATGGATAGAAGGACTTTCA GTGTTTCTGCGAGGGACCTTGgatgaaaaaattaaatATTGTTTTGAGGTTTATGACTTA AATGGGGATGGATATATTTCACGAGAAGAGATGTTCCACATGCTGAAAAACAGTCTCATA AAACAACCAACAGAAGAAGATCCAGATGAAGGGATTAAGGACTTGGTAGAGATTACTCTT AAAAAGATGgaCCACGATCACGACAGCAGACTTTCATACGCTGATTTTGAGAAAGCAGTA AAAGAAGAAAATCTCTTGCTTGAGGCTTTTGGAGCTTGTCTTCCTGATGCAAAGagtaTT CTTGCTTTTGAAAGACAGGCCTTCCAG------GATACCACAGAAAAT
atgctgaaaatgTCGGCGATGAACAGAAAATTAATTCAAAACCTCGCCGAGACTTTATGC AGACAAGTCAAACATTTTAATAAAACAGAGACGGAGTGTCTGATAAGGCTGTTCAACAGT CTGCTGGGAGAGCAGGCAGAGAGAAAGACGACTATTGGAGTGGACCGGGCCAAATTCAGA AATATACTGCACCACACTTTCGGGATGACCGACGACATGATGACGGACAGAGTTTGTCGT GTCATTGACAAGGACAACGATGGCTACTTAAGCGTTAAAGAGTGGGTTGAggctctgtct gtctttctaagAGGCACACTGGATGAAAAAATGAAATaCTGTTTTGAGGTGTATGACCTG AACGGGGATGGATACATCTCACGTGAGGAGATGTTTCAGATGCTGAAAGACAGCCTCATC AGGCAGCCCACCGAAGAGGATCCTGATGAGGGGATTAAGGATATTGTGGAGATTGCCTTG AAAAAAATGGATTATGACCATGATGGAAGAGTTTCTTATGCTGATTTTGAGAAGACGGTC ATGGATGAAAACCTTTTACTAGAAGCTTTTGGAAACTGCCTTCCTGATGCAAAGAGTGTA CTAGCATTTGAGCAACAGGCATTCCAGAAACACGAACACTGCAAAGAA
seq3-seq4
ATGGATCGCCATTCCAAtttaatttccatttggctgcaACTGGAACTGTGTGCCATGGCA GTACTTCTGGCAAAAGGGGAGATAAGATGCTACTGTGATGCAGCGCATTGTGTGGCAACA GGTTACATGTGTAAATCCGAGTTAAATGCCTGCTTCACCAGGCTTCTGGACCCACAGAAC ACAAACTCCCCTCTCACGCATGGCTGCTTGGACCCGACTGCAAACACAGCAGATGTTTGC CATGCTGGAAGGACAGAGAGCCGCGCTGGGGCCTCGGAGAAGCTTGAGTGCTGTCACGAC GATATGTGCAATTACAGAGGACTCCATGATGTTGTTTCATATCCCAGGGGGGACAGCTCA GATCATGGAACAAGATATCAGCCAGACAGTAGCAGGAATCTTCTGACCAGGGTTCAGGAT TTAACATCCTCTAAAGAGCTGTGGTTCAGAGCAGCCGTGATCGCTGTGCCCATCGCTGGG GGGCTCATTCTAGTGCTTCTCATCATGCTCGCCTTGCGGATGCTTCGAAGTGAAAACAAA AGACTGCAGGACCAGAGGCAGCAGATGCTGTCCCGCTTGCACTACAACTTTCATGGA--CACCACACGAAGAAGGGCCAGGTAGCCAAACTGGATTTGGAATGCATGGTTCCCGTAACC GGACACGAGAACTGCTGTATGACTTGCGACAAACTGCGACAGTCTGAACTCCACAAT-----------------GATAAATTGCTGTCTTTAGTTCACTGGGGAATTTACAGCGGTCAC GGGAAATTGGAATttgta
ATGGATCGC---------CTGGTTTCTCTGTGGTTTCAGCTGGAACTTTGTGCGATGGCT GTTCTTCTCACGAAAGGAGAGATCAGGTGCTACTGTGACGCACCGCACTGCGTTGCCACC GGATACATGTGTAAATCAGAGCTCAACGCTTGCTTTACTAAGGTCCTGGACCCTCTTAAC ACAAACTCACCTTTAACACACGGCTGCGTGGATTCGCTTTTAAACTCTGCAGACGTGTGC TCTAGTAAAAATGTGGACATTTCAAGTGGAAGCTCCTCTCCTGTGGAGTGCTGCCATGAT GATATGTGTAACTACAGGGGTTTGCATGAC---CTCACACACCCCAGAGGGGACTCAACA GAC---------CGATACCACAGC---TCCAATCAGAACCTGATCACAAGGGTGCAAGAG TTAGCGTCTGCTAAAGAGGTGTGGTTCCGGGCGGCGGTGATAGCGGTTCCCATCGCGGGT GGGCTTATCCTGGTTCTGCTGATTATGCTGGCGTTGCGAATGCTCCGTAGCGAAAACAAG CGTCTCCAGGCACAGCGCCAGCAGATGCTTTCTCGCCTGCATTACAGCTTTCACGGACAC CACCATGCCAAGAAAGGCCACGTGGCTAAGTTGGACTTGGAGTGTATGGTGCCGGTAACG GGACATGAGAACTGTTGTCTGGGCTGCGATAAGCTGCGGCAGACGGATTTGTGCACTGGA GGAGGAAGCGGGGGTGAGCGTCTCCTATCTCTGGTACACTGGGGGATGTACACGGGGCAC GGAAAGCTGGAGTTCGTA
```

## Batch 4DTV calculation
```
batch_4DTV_calculation.pl Merged.AXT > Merged.4DTV
```
The 4DTV results are in Merged.4DTV

