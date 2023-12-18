##Downloaded 30. Aug. 2022

# archea virus fungus
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/*genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/*genomic.fna.gz

gzip -d *
cat *.fna > archaea_fungi_virus.fa
rm *.fna
bowtie2-build --threads 50 archaea_fungi_virus.fa archaea_fungi_virus.fa

#vert_other
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_other/*genomic.fna.gz
gzip -d *
cat *.fna > vert_other.fa
rm *.fna
perl /home/tvg137/fasta-splitter.pl --n-parts 8 vert_other.fa
rm vert_other.fa
i=1
for file in vert_other.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in vert_other.?
do
bowtie2-build --threads 50 $file $file
done

#vert_mam
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/*genomic.fna.gz
gzip -d *
cat *.fna > vert_mam.fa
rm -f *.fna
perl /home/tvg137/fasta-splitter.pl --n-parts 10 vert_mam.fa
rm vert_mam.fa

i=1
for file in vert_mam.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in vert_mam.*
do
bowtie2-build --threads 50 $file $file
done


# invert
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/invertebrate/*genomic.fna.gz
gzip -d *
cat *.fna > invert.fa
rm *.fna
perl /home/tvg137/fasta-splitter.pl --n-parts 3 invert.fa
rm invert.fa
i=1
for file in invert.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in invert.?
do
bowtie2-build --threads 50 $file $file
done

# plant
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/*genomic.fna.gz

gzip -d *
cat *.fna > plant.fa
rm *.fna
perl /home/tvg137/fasta-splitter.pl --n-parts 5 plant.fa
rm plant.fa

i=1
for file in plant.part*
do
bname=$(basename "$file" | cut -d. -f1)
mv $file $bname.$i
i=$(expr ${i} + 1)
done
for file in plant.?
do
bowtie2-build --threads 50 $file $file
done

# mito
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/*genomic.fna.gz #1
gzip -d *
cat *.fna > mitochondrion.fa
rm *.fna
bowtie2-build --threads 50 mitochondrion.fa mitochondrion.fa

# protozoa
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protozoa/*genomic.fna.gz #1
gzip -d *
cat *.fna > protozoa.fa
rm *.fna
bowtie2-build --threads 50 protozoa.fa protozoa.fa


# plastid
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/*genomic.fna.gz 
gzip -d *
cat *.fna > plastid.fa
rm *.fna
bowtie2-build --threads 50 plastid.fa plastid.fa
