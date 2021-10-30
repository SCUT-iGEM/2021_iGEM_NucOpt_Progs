function fastaWrite(pop,path)

N=length(pop);
fastaFormat={};
for i=1:N
    str=pop{i};
    len=size(str,2);
    fastaStr="<sequence ";
    for j=1:70:len
        last=min(j+69,len);
        fastaStr=sprintf(fastaStr+newline+str(j:last));
    end
    fastaFormat=[fastaFormat;fastaStr];
end
xlswrite(path, fastaFormat, 'Sheet1', 'A1')
end

