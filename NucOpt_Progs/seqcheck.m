function isok=seqcheck(seq,parent,runsstart)%使用全局变量减少参数传递
%The sole purpose of this program is to make sure that a sequence can be
%synthesized by IDT's gblocks.  It was sufficient at the time of writing but some
%features of it may no longer be necessary as synthesis technology
%improves.

%initiates things
isok=1;
complement=seqcomplement(seq);
len=length(seq);
GCstart=gcprofile(parent);
forbiddenruns={'AAAAAAAAA','CCCCCC','TTTTTTTTT','GGGGGG'};

%check total GC content
GC=gccontent(seq);
if GC>75||GC<=25
    isok=0;
end

%check GC content every 100bp
if isok==1
    GC=gcprofile(seq);
    if max(GC)>80||min(GC)<24 %checks if GC content is not within acceptable range
        if min(GC)<min(GCstart) %is minimum GC content of new sequence lower than parent?
            isok=0;
        elseif min(GC)==min(GCstart) %is minimum GC content of new sequence equal to parent?
%%%更改     if sum(GC==min(GC))>=sum(GCstart==min(GCstart)) %is there not less of the minimum GC value than for the parent?
            if sum(GC==min(GC))>sum(GCstart==min(GCstart))
                isok=0;
            end
        end
    end
end

%Check for homopolymers which are too long
if isok==1    
%%%更改%%
    seqforbidden=containsforbidden(seq,forbiddenruns);
    if size(cell2mat(seqforbidden),2)>size(cell2mat(runsstart),2)
       isok=0; 
    end
%
%     seqforbidden=containsforbidden(seq,forbiddenruns);
%     if isempty(seqforbidden{1})==0||isempty(seqforbidden{2})==0||isempty(seqforbidden{3})==0||isempty(seqforbidden{4})==0 %Are there runs?
%         for n=1:4
%             if isempty(seqforbidden{n})==1
%                 moreruns(n)=0;
%             else
%                 moreruns(n)=(sum(seqforbidden{n})<sum(runsstart{n}));%非空则
%             end
%         end
%         if sum(moreruns)==0 % Are there equal or more runs in the new sequence than the parent sequence?
%                             % equal不行啊？？？
%             isok=0;
%         end
%     end
end
            

% check for hairpins
if isok==1
    %create the dot matrix and rotate it
    for n=1:len
        hpdot(n,:)=complement==seq(n);
    end
    hpdot=rot90(hpdot);
    %search through all the diagonals for runs of "dots" of a certain
    %length.  higher than 8bp hairpin is bad if GC content is greater than
    %80%, higher than 11bp hairpin is always bad
    for n=-length(seq)+1:length(seq)-1
        if isok==1;
            a=diag(hpdot,n);
            dia=a(1:ceil(length(a)/2));
            hpsmall=strfind(dia',ones(1,9));
            hpbig=strfind(dia',ones(1,12));
            if size(hpbig)~=0
                isok=0;
            elseif size(hpsmall)~=0
                for m=1:size(hpsmall)
                    if n<=0
                        GC=gccontent(seq(hpsmall(m):hpsmall(m)+6));
                    else
                        GC=gccontent(seq(len-hpsmall(m)-5:len-hpsmall(m)+1));
                    end
                    if GC>80
                        isok=0;
                    end
                end
            end
        end
    end
end

%check for repeats or inverted repeats in each 100bp subsequence
% done=0;
% if isok==1
%     %create the dot matrices
%     for m=1:len
%         repdot(m,:)=seq==seq(m);
%     end
%     %straighten out the diagonals for only those nucleotides within a 100bp
%     %window of one another
%     for m=1:95
%         repdiags(:,m)=padarray(diag(repdot,m),length(diag(repdot,1))-length(diag(repdot,m)),'post');
%     end
%     %iterate through all relevant repeat lengths
%     for l=4:50
%         reps=zeros(100-l,len);
%         if isok==1&&done==0
%             %look for the relevant repeat length within a 100bp window and
%             %save all instances of the repeat to reps as its location (if
%             %the repeat would extend beyond a 100bp window then it isn't
%             %counted, hence the 100-l)
%             for m=1:100-2*l
%                 occur=strfind(repdiags(:,m+l-1)',ones(1,l));
%                 if size(occur,1)>0
%                     reps(m,:)=padarray(occur,[0 len-length(occur)],'post');
%                 end
%             end
%             %and filter out any repeats which occur within l bp of one
%             %another
%             if sum(sum(reps))>0
%                 filteredreps=reps;
%                 for n=1:100-2*l+1
%                     a=~ismember(filteredreps(n+1:n+l-1,:),filteredreps(n,:));
%                     filteredreps(n+1:n+l-1,:)=a.*filteredreps(n+1:n+l-1,:);
%                 end
%                 res = [(1:len)' histc(filteredreps(:), 1:len)];
%                 sortedres = sortrows(res, -2);
%                 %check if the repeats cover more than 30bp (ie 30% of the 100bp
%                 %window)
%                 if (sortedres(1,2)+1)*l>=30&&sortedres(1,2)>0
%                     isok=0;
%                 end
%             else
%                 done=1;
%                 %stops calculation early if it can be.
%             end
%         end
%     end
% end