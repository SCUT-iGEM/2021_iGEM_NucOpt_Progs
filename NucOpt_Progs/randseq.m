function sequence=randseq(n,gc)
%makes a random sequence of the length specified by n and with an average
%gc content of gc
randnum=randi(100,[1,n]);
for i=1:n
    if randnum(i)>=1&&randnum(i)<gc/2
        sequence(i)='C';
    elseif randnum(i)>=gc/2&&randnum(i)<gc
        sequence(i)='G';
    elseif randnum(i)>=gc&&randnum(i)<(gc+(100-gc)/2)
        sequence(i)='T';
    else
        sequence(i)='A';
    end
end