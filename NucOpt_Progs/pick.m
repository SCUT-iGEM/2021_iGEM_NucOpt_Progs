function tomutate = pick(seqnum,numchanges)
length=size(seqnum,2)-numchanges+1;
tomutate=zeros(length,numchanges);
for j=1:numchanges
    tomutate(1,j)=seqnum(1)+j-1;
end

for i=2:length-numchanges+1
    tomutate(i,:)=tomutate(i-1,:)+1;
end
end

