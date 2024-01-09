# Randseq
A tool to generate toy random sequences.
For instance, these command can be used to generate toy sequences for kreeq:
```
./build/bin/randseq --verbose -l 10000 --gc-content 0.2 --mutation-rate 0.01 --average-read-length 50 --read-coverage 50
kreeq validate -f referenceError.fasta -r reads.fastq -o test.bed
awk 'NR==FNR{a[$1$2]=1;next}{if (a[$1$2] == 1) {b=0} else {b=1}; printf $0"\t"b"\n";}' errors.bed test.bed > test_10000.bed
```
