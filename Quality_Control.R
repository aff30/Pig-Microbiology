### FastQC quality report
#This process was done in ROAR - Penn State's supercomputer

# Activate conda environment
conda activate bioinfo

# make a directory to store the fastQC report
mkdir quality

# check your files
ls

# Run fastqc
fastqc *.fastq.gz -o quality

# If you encounter an error when running multiQC, you should run the following commands:
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

#Run MultiQC - this will generate a html report
multiqc quality/* -o .

# if you are in a super computer, you can transfer your html folder to your local computer using the scp function (google scp and you will find a way on how to transfer because it place is one case)
#once it is local you can open the file
#Example:
scp -r credentials@submit.aci.ics.psu.edu:/gpfs/group/evk5387/default/pig_usp/multiqc_data .
scp -r credentials@submit.aci.ics.psu.edu:/gpfs/group/evk5387/default/pig_usp/multiqc_report.html .
