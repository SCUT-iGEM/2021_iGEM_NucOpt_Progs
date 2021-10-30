function [maxsequence,maxarea]=nucleomin(sequence,prombeg,promend,numchanges,forbiddensites,forbiddenseqs)
%nucleomin takes a sequence as input and computes the n-nucleotide variant
%with the minimum nucleosome affinity that is also synthesizable and also 
%does not have any additional or fewer transcription factor binding sites. n is user-defined.  

%input sequence must be uppercase strings
%ï¿½ï¿½ï¿½ï¿½ï¿½seqï¿½ï¿½ï¿½ï¿½ï¿½Ç´ï¿½Ð´ï¿½ï¿½Ä¸

%forbiddenseqs must be a cell array with each motif in the first row.
%motifs specified in forbiddenseqs will neither be created or destroyed.
%This is for things like TATA boxes or other general purpose transcription
%factors which may be present.  Also ATGs if you like.

%forbiddensites is a row vector of positions that you don't want the
%program to mutate.  For example, things like transcription factor binding
%sites.

%numchanges tells the program how far to search from the parent sequence to
%find an improved promoter.  numchanges=1 searches all single mutants,
%numchanges=2 searches all double mutants, etc... 

%Prombeg and promend specify the positions of the beginning and end of the
%promoter in "sequence"  We recommend prombeg be at least 200.

forbiddenruns={'AAAAAAAAA','CCCCCC','TTTTTTTTT','GGGGGG'};
runsstart=containsforbidden(sequence(prombeg:promend),forbiddenruns);
%IDT doesn't like these sequences, so we're making note of where they are.
%For determining if a sequence is synthesizable

maxaffinities=affinity(sequence);
learningcurve=maxaffinities(1:25);
%This is the affinity of the sequence we're starting from.

refforbidden=containsforbidden(sequence,forbiddenseqs);%ï¿½ï¿½Ê¼Ê±ï¿½ï¿½Ö¹Í»ï¿½ï¿½ï¿½Î»ï¿½ï¿?
%We also save the locations of anything in forbiddenseqs.  For determining
%if a sequence contains any extra transcription factor binding sites.

%%tomutate=pick(prombeg:promend,numchanges,'r');%tomutateï¿½ï¿½Â¼ÒªÍ»ï¿½ï¿½ï¿½Î»ï¿½Ã£ï¿½ï¿½ï¿?
%%%pickï¿½ï¿½ï¿½ï¿½ï¿½Ä¹ï¿½ï¿½ï¿½ï¿½ï¿½Ê²Ã´ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ê²Ã´ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
% % ï¿½Ô¼ï¿½Ð´ï¿½ï¿½pickï¿½ï¿½
tomutate=pick(prombeg:promend,numchanges);

%%ï¿½ï¿½ï¿½ï¿½%%
% % tomutate=[prombeg:promend]';


%%ï¿½ï¿½ï¿½ï¿½%%basechanges=str2digit(dec2base(0:4^numchanges-1,4,numchanges));%ï¿½ï¿½ï¿½ï¿½4ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ê¾ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½nÎ»
basechanges=str2num(dec2base(0:4^numchanges-1,4,numchanges));
%generates a worklist for all the bases to mutate during the search
%for an improved promoter.  For each entry in tomutate, basechanges is a
%worklist for what to mutate those bases to in its search.

n=1;
testarea=[];
tic
%n and tic are just there if you are impatient and want to see the progress
%of nucleomin in real time. Also initializing testarea.

for i=1:size(tomutate,1)
    %cycles through all the positions needing to be randomized
    
    for j=1:size(basechanges,1)
        %cycles through all possible bases at the randomized positions
        
        badseq=0;
        %badseq is 1 if sequence has an issue and should be thrown out,
        %badseq is 0 otherwise.
        
        testseq=sequence;
        %this is the sequence we're going to be mutating
        
        unicom=[tomutate(i,:)' basechanges(j,:)'];
        
        if length(unique(tomutate(i,:)'))==size(unique(unicom,'rows'),1)%È·±£bpÊý´óÓÚ1Ê±Ã»ÓÐ²»±ØÒªµÄÖØ¸´
            % the previous two lines are for making sure that for more than
            % nbp mutations at a time (n>1), that the (<n)bp mutants are
            % also computed and without unnecessary repetitions
            
            for k=1:numchanges
                %making the specified mutations to testseq
                
                
                if sum(cell2mat(forbiddensites)==tomutate(i,k))==0 %È·±£½ûÖ¹Í»±äµÄÐòÁÐ²»Í»±ä
                %%¸ü¸Ä%% 
                %%if okMutate(forbiddensites,forbiddenseqs,tomutate(i,:))
                %makes sure we're not going to mutate anything in
                    %forbiddensites
                    
                    if basechanges(j,k)==0
                        if sequence(tomutate(i,k))=='A'
                            badseq=1;
                            %prevents us from mutating to the same base, which
                            %would eat up time.
                        else
                            testseq(tomutate(i,k))='A';
                            %make the mutation
                        end
                    elseif basechanges(j,k)==1
                        if sequence(tomutate(i,k))=='C'
                            badseq=1;
                        else
                            testseq(tomutate(i,k))='C';
                        end
                    elseif basechanges(j,k)==2
                        if sequence(tomutate(i,k))=='T'
                            badseq=1;
                        else
                            testseq(tomutate(i,k))='T';
                        end
                    elseif basechanges(j,k)==3
                        if sequence(tomutate(i,k))=='G'
                            badseq=1;
                        else
                            testseq(tomutate(i,k))='G';
                        end
                    end
                else
                    badseq=1;
                end
            end
        else
            badseq=1;
        end
        if badseq==0
            testforbidden=containsforbidden(testseq,forbiddenseqs);%ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ö¹Í»ï¿½ï¿½ï¿½Î»ï¿½ï¿½
            %looks for forbidden motifs in the mutated sequence
            
            isok=seqcheck(testseq(prombeg:promend),sequence(prombeg:promend),runsstart);
            %makes sure IDT can synthesize the mutated sequence.
            
            if isok==1
% %                 pause();
            end
            
            
            if isequal(refforbidden,testforbidden)==0||isok==0
                badseq=1;
                %A sequence is bad if it contains a different number of
                %forbidden motifs than the starting sequence or if it
                %cannot be synthsized by IDT.
            end
        end
        if badseq==0;
            testaffinity=affinity(testseq);
            testarea(n)=seqarea(learningcurve,testaffinity,prombeg-73,promend-73);         
            testseqs{n}=testseq;
            %if the sequence is ok, this will compute the nucleosome
            %affinity under the promoter and add this area to the list of
            %mutants
            n=n+1;
           percentdone=((i-1)*size(basechanges,1)+j)/(size(tomutate,1)*size(basechanges,1))*100
           timeleft=toc/(percentdone/100)-toc
%   you can enable the previous two lines if you are impatient and want to
%   see progress of nucleomin.
        end
    end
end
[maxarea,i]=min(testarea);%
maxsequence=testseqs{i};
%finds the promoter with the minimum nucleosome affinity and returns it.
toc