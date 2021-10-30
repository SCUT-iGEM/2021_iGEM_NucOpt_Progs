function sequencefix=remforbidden(sequence,forbiddensites,forbiddenseqs)
%Tries to remove as many motifs found in forbiddenseqs as possible from
%sequence given that no bases in forbiddensites can be changed

isbetter=1;
tic
while isbetter==1;
    problemsites=problemrank(sequence,forbiddenseqs);
    %notes which bases contain motifs in forbiddenseqs
    isbetter=0;
    for i=problemsites
        %iterates through the problem bases
        for j=0:3
            %iterates through all basepair changes
            if isbetter==0;
                % if we haven't removed a motif yet
                badseq=0;
                testseq=sequence;
                if sum(forbiddensites==i)==0
                    %makes sure we aren't mutating a forbidden site
                    if j==0
                        if testseq(i)=='A'
                            badseq=1;
                        else
                            testseq(i)='A';
                        end
                    elseif j==1
                        if testseq(i)=='C'
                            badseq=1;
                        else
                            testseq(i)='C';
                        end
                    elseif j==2
                        if testseq(i)=='T'
                            badseq=1;
                        else
                            testseq(i)='T';
                        end
                    elseif j==3
                        if testseq(i)=='G'
                            badseq=1;
                        else
                            testseq(i)='G';
                        end
                    end
                else
                    badseq=1;
                end
                if badseq==0
                    testsites=problemrank(testseq,forbiddenseqs);
                    %counts the forbidden motifs in the new sequence
                    if length(testsites)<length(problemsites)
                        isbetter=1;
                        sequence=testseq;
                        %if the new sequence contains less motifs than the
                        %original, discard the parent and save the good one
                        %for the next round
                    end
                end
            end
        end
    end
end
sequencefix=sequence;
%this sequence should contain the minimum number of forbidden motifs given
%that we can't change anything in forbiddensites.
toc