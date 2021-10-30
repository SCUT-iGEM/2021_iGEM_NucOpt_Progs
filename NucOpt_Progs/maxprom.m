function [proms,areas,curves]=maxprom(sequence,prombeg,promend,numchanges,forbiddenseqs)
%This program will take a promoter and iteratively make improvements to it
%in 'numchanges' increments until it can no longer be improved.  See the
%functions this program calls for what each of the input arguments are.

%Outputs the promoter, summed affinity area, and the affinity values along
%the promoter for each iteration in the optimization

global RPath;
global dirPath;
RPath='G:\R\R\R-4.1.0\bin\R.exe';
dirPath='G:\iGEM\proNuGO';

% -----------继续跑取消注释下面一段---------
% load('maxpromdata.mat');
% n=size(proms,2);
% -----------------------------------------
 

% % -----------重新跑取消注释下面一段---------
t_loop=[];%记录每个循环耗费的时间
rewriteR_NuPoP(promend);%重写预测用的.R文件
forbiddensites=containsforbidden(sequence,forbiddenseqs);
proms{1}=sequence;
refaffinity=affinity(sequence);
curves(1,:)=refaffinity;
learningcurve=refaffinity(1:25);
areas(1)=seqarea(learningcurve,refaffinity,prombeg-73,promend-73);
save('maxpromdata.mat')
%initializes things
t1=clock
[proms{2},~]=nucleomin(sequence,prombeg,promend,numchanges,forbiddensites,forbiddenseqs);
testaffinity=affinity(proms{2});
[areas(2),curves(2,:)]=seqarea(learningcurve,testaffinity,prombeg-73,promend-73);

t2=clock
t_loop=[t_loop;etime(t2,t1)];
save('maxpromdata.mat')
%does the first iteration of promoter improvement

t1=clock;
[proms{3},~]=nucleomin(proms{2},prombeg,promend,numchanges,forbiddensites,forbiddenseqs);
testaffinity=affinity(proms{3});
[areas(3),curves(3,:)]=seqarea(learningcurve,testaffinity,prombeg-73,promend-73);

t2=clock;
t_loop=[t_loop;etime(t2,t1)];
save('maxpromdata.mat')
%and the second

t1=clock;
[proms{4},~]=nucleomin(proms{3},prombeg,promend,numchanges,forbiddensites,forbiddenseqs);
testaffinity=affinity(proms{4});
[areas(4),curves(4,:)]=seqarea(learningcurve,testaffinity,prombeg-73,promend-73);

t2=clock;
t_loop=[t_loop;etime(t2,t1)];
save('maxpromdata.mat')
%and the t hird

n=4
% -----------------------------------------


while isequal(proms{n},proms{n-2})==0
    %continues to make improvements while we aren't stuck in a local
    %minimum.
    n
    n=n+1;
    t1=clock;
    [proms{n},~]=nucleomin(proms{n-1},prombeg,promend,numchanges,forbiddensites,forbiddenseqs);
    testaffinity=affinity(proms{n});
    [areas(n),curves(n,:)]=seqarea(learningcurve,testaffinity,prombeg-73,promend-73);
    t2=clock;
    t_loop=[t_loop;etime(t2,t1)];
    %makes improvements
    
    areas
    %things without semicolons are just so you can see how far along it is.
    % Feel free to take these out if you're hardcore.
    
    save('maxpromdata.mat')
    % just in case your computer crashes you can start from where you left
    % off.
end
etime(t2,t1)
