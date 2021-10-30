function [sequence,forbiddensites]=randprom(params,gc)
%makes an intial random sequence for the synthetic promoter and generates
%forbiddensites from the nonrandom (user-specified) regions of params.

forbiddensites=[];
sequence=[];
for n=1:length(params)
    if ischar(params{n})==1
        a=length(sequence)+1;
        sequence=strcat(sequence,params{n});
        b=length(sequence);
        forbiddensites=[forbiddensites a:b];
    else
        sequence=strcat(sequence,randseq(params{n},gc));
    end
end