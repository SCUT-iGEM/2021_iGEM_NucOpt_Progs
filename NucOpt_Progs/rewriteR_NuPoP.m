function rewriteR_NuPoP(promend)
global dirPath;
text = fileread("../../R/NuPoPTemplate.txt");
text=sprintf(text,promend+101-73-73,promend+101,promend+101-73);%outputά�ȣ�Ԥ��ά�ȣ�����ά��
Rpath=replace([dirPath,'\IOfile'],'\','/');
text=sprintf(text,Rpath);
fid=fopen("../../R/NuPoP.R",'wb');
fwrite(fid,text);
fclose(fid);
end

