function [area,affinitiesscaled]=seqarea_0(learningcurve,affinities,prombeg,promend)
%scales an input affinity curve so that it is comparable to another
%affinity curve, specified by 'learningcurve'

%Npred does this weird thing where it scales its affinity output depending
%on the max and min affinities.  This becomes a problem when promoter
%modifications change the extrema of the affinity curve.  This program
%effectively "unscales" this affinity curve and "rescales" it to the same
%scale as a reference curve.  To do so, it does a simple linear regression
%between the first 25 affinity values of each curve.  In order for this to
%be valid, the first 25bp of each generating seqence must be the same.  We
%recommend putting 200bp of context before and 100bp of context after each
%promoter, as detailed in "Notes for MATLAB".  prombeg and promend specify
%the locations of the first and last nucleotide of the promoter of interest
%in 'affinities', for calculation of areas.  

%%？？
%Npred做了一件奇怪的事情，它根据最大和最小亲和度来缩放亲和度输出。当启动子
%修饰改变亲和曲线的极值时，这就成了一个问题。这个程序有效地“取消”这个亲
%和曲线，并将它“重新缩放”到与参考曲线相同的比例。为此，它在每条曲线的前
%25个亲和值之间进行简单的线性回归。为了使其有效，2每个生成序列的前5bp必须相同。
%我们建议在每个启动子之前放置200bp的上下文，在每个启动子之后放置100bp的上下文，
%详见“MATLAB注释”。prombeg和promend指定“亲和力”中感兴趣的启动子的第一个和
%最后一个核苷酸的位置，用于计算面积。



[~,m,b]=regression(affinities(1:25)',learningcurve');
%Performs regression between learningcurve and affinities

affinitiesscaled=affinities.*m+b;
%Uses regression to scale affinities


area=sum(affinitiesscaled(prombeg:promend));
%computes nucleosome area under promoter region.