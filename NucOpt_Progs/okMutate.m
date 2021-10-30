function flag = okMutate(forbiddensites,forbiddenseqs,tomutate)
flag=1;
forbiddensizes=zeros(size(forbiddenseqs));
for i=1:size(forbiddenseqs,2)
    forbiddensizes(i)=size(forbiddenseqs(i),2);
end
for i=1:size(forbiddensites,2)
    for j=1:size(forbiddensites{i},2)
        if((tomutate(1)>=forbiddensites{i}(j))&&(tomutate(1)<=forbiddensites{i}(j)+forbiddensize(i))...
            &&(tomutate(end)>=forbiddensites{i}(j))&&(tomutate(end)<=forbiddensites{i}(j)+forbiddensize(i)))
            flag=0;
            return;
        end
    end
end
return;
end

