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

%%����
%Npred����һ����ֵ����飬������������С�׺Ͷ��������׺Ͷ��������������
%���θı��׺����ߵļ�ֵʱ����ͳ���һ�����⡣���������Ч�ء�ȡ���������
%�����ߣ����������������š�����ο�������ͬ�ı�����Ϊ�ˣ�����ÿ�����ߵ�ǰ
%25���׺�ֵ֮����м򵥵����Իع顣Ϊ��ʹ����Ч��2ÿ���������е�ǰ5bp������ͬ��
%���ǽ�����ÿ��������֮ǰ����200bp�������ģ���ÿ��������֮�����100bp�������ģ�
%�����MATLABע�͡���prombeg��promendָ�����׺������и���Ȥ�������ӵĵ�һ����
%���һ���������λ�ã����ڼ��������



[~,m,b]=regression(affinities(1:25)',learningcurve');
%Performs regression between learningcurve and affinities

affinitiesscaled=affinities.*m+b;
%Uses regression to scale affinities


area=sum(affinitiesscaled(prombeg:promend));
%computes nucleosome area under promoter region.