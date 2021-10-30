function [proms,areas,curves]=synthprom(params,prombeg,promend,numchanges,forbiddenseqs)
%this function takes a general outline for a promoter (specified in params)
%and makes a synthetic nucleosome optimized promoter.  All variables are
%the same for nucleomin except for params.  Params is a cell array whose
%contents are either numbers or DNA sequences.  numbers represent length of
%random (nucleosome optimizable) unspecified DNA sequences, and DNA
%sequences are TFBSs or anything else you want to keep constant during the
%optimization.  Put each segement in the order you want it to appear.

%Example: {{3},{'AGTAGCA'},{7}} is NNNAGTAGCANNNNNNN

gc=35;
%GC content of randomly generated portions of the promoter.  Yeast is
%around 35% but if you are in a different organism you can change that
%here.

[sequence,forbiddensites]=randprom(params,gc);
%makes an intial random sequence for the synthetic promoter and generates
%forbiddensites

sequencefix=remforbidden(sequence,forbiddensites,forbiddenseqs);
%removes anything in forbiddenseqs (like TFBSs) randomly generated in sequence, unless
%they're contained in forbiddensites.

[proms,areas,curves]=maxprom(sequencefix,prombeg,promend,numchanges,forbiddensites,forbiddenseqs);
%Takes the initialized sequence and performs a nucleosome optimiztion.