
srun --export=ALL --ntasks-per-node 24  --nodes 1 --partition=comppriority --mem 65G --pty bash -i commands1.sh

srun --export=ALL --ntasks-per-node 5  --nodes 1 --partition=comptest --mem 5G --pty bash


# commands generator

#!/bin/bash

# Define the input files
sample_list="sample.list"
database_list="database.list"

# Loop through sample names
while IFS= read -r file; do
    # Loop through database names
    while IFS= read -r DB; do
        # Generate the command
        command="bowtie2 -k 1000 -x $DB -U $file --no-unal --threads 24  -S $file.$(basename $DB).sam"
        # Print the command to a text file
        echo "$command" >> commands.txt
    done < "$database_list"
done < "$sample_list"




##### slurm command bash script (sbatch script) OBS change array 1-46 should be the number of lines in the command.txt file

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=00:10:00
#SBATCH --mem=65G
#SBATCH --array=1-1%100
#SBATCH --export=ALL
#SBATCH --output=holiMap_%A_%a.out

set -x
source /home/npl206/.bashrc
cd /projects/lundbeck/people/npl206/brick_aDNA
conda activate KapKBH
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
task=($(awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}' /projects/lundbeck/people/npl206/brick_aDNA/commands.txt))
"${task[@]}"


### execute the sbatch file

sbatch slurm_executerHoli.sh


### database.list file
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.1
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.2
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.3
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.4
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.5
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.6
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.7
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.8
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.9
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.10
/projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.1
/projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.2
/projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.3
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.1
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.2
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.3
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.4
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.5
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.6
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.7
/projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.8
/projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.1
/projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.2
/projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.3
/projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.4
/projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.5
/projects/wintherpedersen/data/refseq_30Aug2022/archaea_fungi_virus/archaea_fungi_virus.fa
/projects/wintherpedersen/data/refseq_30Aug2022/plastid/plastid.fa
/projects/wintherpedersen/data/refseq_30Aug2022/mitochondrion/mitochondrion.fa
/projects/wintherpedersen/data/refseq_30Aug2022/protozoa/protozoa.fa
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.1
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.2
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.3
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.4
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.5
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.6
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.7
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.8
/projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.9
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.1
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.2
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.3
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.4
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.5
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.6
/projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.7


31G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.1
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.2
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.3
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.4
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.5
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.6
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.7
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.8
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.9
32G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_mammal/vert_mam.10
35G /projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.1
36G /projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.2
35G /projects/wintherpedersen/data/refseq_30Aug2022/invertebrate/invert.3
39G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.1
39G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.2
39G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.3
39G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.4
38G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.5
37G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.6
38G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.7
37G /projects/wintherpedersen/data/refseq_30Aug2022/vertebrate_other/vert_other.8
22G /projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.1
23G /projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.2
24G /projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.3
23G /projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.4
22G /projects/wintherpedersen/data/refseq_30Aug2022/plant/plant.5
9.1G /projects/wintherpedersen/data/refseq_30Aug2022/archaea_fungi_virus/archaea_fungi_virus.fa
716M /projects/wintherpedersen/data/refseq_30Aug2022/plastid/plastid.fa
187M /projects/wintherpedersen/data/refseq_30Aug2022/mitochondrion/mitochondrion.fa
1.8G /projects/wintherpedersen/data/refseq_30Aug2022/protozoa/protozoa.fa
44G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.1
45G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.2
45G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.3
44G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.4
45G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.5
45G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.6
45G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.7
44G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.8
45G /projects/wintherpedersen/data/ncbi_nt_Sept2022/nt.9
9.8G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.1
9.9G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.2
9.9G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.3
9.9G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.4
9.9G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.5
9.9G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.6
9.9G /projects/wintherpedersen/data/norway_plant_ctgenoms/complete_genomes/final_db/norPlantCom.7
