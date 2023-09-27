##Reverse reads had bad quality. We tried running fastp (a tool that trims the reads and remove NAs) with them and even after fastp DADA2 did not work! Sequences were not merging for some reason!
# For more information on fastp, access: https://github.com/OpenGene/fastp
#Here is a summary of what was done:

# First create a list with the name of your samples
# you can create this list in the visual studio code and then same in your directory

# 1- copy and past your samples name in the visual code
# remove the _L001_R1_001.fastq.gz and _L001_R2_001.fastq.gz

# 2- remove ao the duplicates (since I have the same names for forward and reverse)

Control + F.
# Toggle “Replace mode”
# Toggle “Use Regular Expression” (the icon with the .* symbol)
# In the search field, type:
^(.*)(\n\1)+$
  # In the “replace with” field, type:
  $1.
#Click. (“Replace All”).

# This will give you a list like the one below: (this list below is an example on how your list should be)
1-base_S105
100_S206
101_S185
102_S100
103_S214
104_S60
105_S46
106_S53
107_S180
108_S103

# Without the files name: _L001_R1_001.fastq.gz and _L001_R2_001.fastq.gz
##YOU CAN DO THIS ON EXCEL AS WELL

# create a new file for the files names
touch filesname
nano filesname
copy and paste the list names you just create

# run fastp using parallel
# if you don't have fastp and parallel you have to download it into to conda environment
# you can check by activating conda environment and check the list

activate conda environment
conda list

# Install newer fastp version
NAMES=filesname
cat $NAMES | parallel /storage/home/aff30/work/fastp -i {}_L001_R1_001.fastq.gz -I {}_L001_R2_001.fastq.gz -o {}_r1_fastp.fq -O {}_r2_fastp.fq -l 150 --cut_window_size 4 --cut_mean_quality 30 --cut_tail --n_base_limit 0
# parallel it is “a shell tool for executing jobs in parallel using one or more compute nodes”

# fastp generates a an html file,transfer it local
scp -r aff30@submit.aci.ics.psu.edu:/gpfs/group/evk5387/default/broiler_microbiome_group/broiler_microbiome_groups/fastp.html


# check reads stat
seqkit stat *_fastp.fq -T -o fastpstat

# read fastp stat file
cat fastpstat