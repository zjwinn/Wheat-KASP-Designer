#!/bin/bash

#SBATCH --job-name="zjwinn"                       #name of the job submitted
#SBATCH -p atlas                                  #name of the queue you are submitting job to
#SBATCH -A gbru_wheat2      
#SBATCH -N 1                                      #number of nodes in this job
#SBATCH -n 48                                     #number of cores/tasks in this job
#SBATCH -t 10:00:00                               #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=zjwinn@ncsu.edu               #email address
#SBATCH --mail-type=BEGIN   
#SBATCH --mail-type=END     
#SBATCH --mail-type=FAIL    
#SBATCH -o "stdout.%x.%j.%N"                      #standard out %x adds job name and %j adds job number to outputfile name and %N adds the node
#SBATCH -e "stderr.%x.%j.%N"                      #optional but it prints our standard error

# Check for local argument
if [ "$1" = "local" ]; then
    # Run test file (local)
    bash ../Wheat-KASP-Designer.sh \
        --input-file 'NC13-20076xGA06493-13LE6_filt_fixed.vcf.gz' \
        --reference-geno  '/mnt/c/Users/zwinn/Music/Ref/161010_Chinese_Spring_v1.0_pseudomolecules.fasta' \
        --snps 'test_snps.txt' \
        --verbose \
        --kasp \
        --debug

else
    # module load miniconda
    module load miniconda3

    # Activate conda environment
    source activate wkd_env

    # Run test file (Atlas)
    bash ../Wheat-KASP-Designer.sh \
        --input-file 'NC13-20076xGA06493-13LE6_filt_fixed.vcf.gz' \
        --reference-geno  '/project/90daydata/gbru_wheat2/zjwinn_project_directory/Ref/161010_Chinese_Spring_v1.0_pseudomolecules.fasta' \
        --snps 'test_snps.txt' \
        --verbose \
        --debug \
        --kasp \
        --max-temp 63 \
        --max-price 200 \
        --max-size 25 
fi




