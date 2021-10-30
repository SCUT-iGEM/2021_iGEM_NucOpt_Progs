function GC=gccontent(seq)
%computes the gc content of seq
numSeq = double(nt2int(seq));
baseNum = [sum(numSeq == 1) sum(numSeq == 2) sum(numSeq == 3) sum(numSeq == 4)];
GC = 100 * ((baseNum(2) + baseNum(3)) / length(numSeq));