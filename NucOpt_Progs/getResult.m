function getResult(Pop,name,prombeg,promend)
N=length(Pop);
%% 画dif表格
allDif=[Pop{1}'];
bpChangeNum=[0];
for i=2:N
    lastref=char(Pop{i-1,1});
    refseq=char(Pop{i,1});
    dif=refseq(:)~=Pop{1}(:);
    bpChangeNum=[bpChangeNum sum(dif)];
    res=[refseq' num2str(dif)];
    allDif=[allDif res];
end

path=['./allDif_',name,'.xlsx'];
xlswrite(path,allDif)
path=['./bpChangeNum_',name,'.xlsx'];
bpChangeNum=[1:N;bpChangeNum];
xlswrite(path,bpChangeNum);
plot(bpChangeNum(1,:),bpChangeNum(2,:));
saveas(gcf,['G:/iGEM/IOfile/Rfile/','bpChangeNum','.jpg']);
%% 画总体图和个体图
mydraw(Pop,prombeg,promend);

end

