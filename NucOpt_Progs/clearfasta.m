function clearfasta(ReadFilePath,WriteFilePath)
%% 2017/06/22 by DQ 
% 此算法为移除文本中空行
FidRead=fopen(ReadFilePath,'rb','ieee-le','UTF-8');
FidWrite=fopen(WriteFilePath,'wb','ieee-le','UTF-8');
while ~feof(FidRead)
    FileRowStr = fgetl(FidRead);
    if ~isempty(FileRowStr )
    fprintf(FidWrite,'%s\n',FileRowStr);
    end
end
fclose(FidRead);
fclose(FidWrite);
delete G:\iGEM\IOfile\sequence1.fasta
end