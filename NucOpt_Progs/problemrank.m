function problemsites=problemrank(sequence,forbiddenseqs)
%Notes the positions of sequence containing a motif found in forbiddenseqs
%and ranks them from lowest nucleotide to highest nucleotide.
badsites=containsforbidden(sequence,forbiddenseqs);
problemsites=[];
for n=1:length(badsites)
    if isempty(badsites{n})==0
        problemsites=[problemsites badsites{n}];
    end
end
problemsites=sort(problemsites);