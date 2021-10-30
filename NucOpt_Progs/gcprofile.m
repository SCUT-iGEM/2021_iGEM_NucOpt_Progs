function GC=gcprofile(seq)
%calculates GC contents in each 100bp sliding window of seq.
len=length(seq);
if len<100
    GC=gccontent(seq);
else
    for n=1:len-99
        GC(n)=gccontent(seq(n:n+99));
    end
end