# DEXSeq_protocol
Create a conda environmanet with python 2.7
Preparing the DEXSeq compatible gff.

    python /home/parashar/anaconda3/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py /home/parashar/archive/rnaseq/gencode.v33.annotation.gtf gencode.v33.gff 2> gff.stderr


DEXSeq doesnot work with BAM files though theyclaim it cam. So convert the bams to sam files using samtools view subcommand.

    while read p
    do
    name=$(basename ${p} .sorted.sam)
    bsub -o $name.o -e $name.e -n 10 -q regularq "python /home/parashar/anaconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -s reverse -a 10 gencode.v33.gff $p ${name}.Readcounts.txt 2> ${name}.stderr"
    done < list

