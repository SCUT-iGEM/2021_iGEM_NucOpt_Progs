# NucOpt Progs Usage

### Please read NuPGO usage tutorial first before reading NucOpt Progs

## Setting

1. Files and path configuration is the same as NuPGO
2. R_path and dir_path setting and its meanings occuring later are the same as both in NuPGO

## Parameters introduction

1. You can find the definition of initial parameters "sequence, prombeg, promend, numchanges, forbiddenseqs" in NuPGO handbook.
2. Notice that numchanges must to be 1, while others will raise error, due to our incomplete repair.

## Operation

1. Initially, please set your Matlab working path as the folder where the 'maxprom' function file exists
2.  Type "max prom(sequence,prombeg,promend,numchanges,forbiddenseqs)" in the command window of Matlab. Of course define those parameter first in the workspace
3. The result will be saved in the working path folder as a file named maxpromdata.mat, whose usage is similar to NuPGO's.
4. t_loop parameter is used for recording the time cost in optimizing each bp.