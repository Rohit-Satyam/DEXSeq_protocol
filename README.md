# DEXSeq_protocol
Create a conda environmanet with python 2.7
Preparing the DEXSeq compatible gff.
We wish to generate nonoverlapping aggregate gtf. For details refer [here](http://seqanswers.com/forums/showthread.php?t=41551). 
> In the process of forming the counting bins, the script might come across overlapping genes. If two genes on the same strand are found
> with an exon of the first gene overlapping with an exon of the second gene, the script’s default behaviour is to combine the genes into a single “aggregate gene” which is subsequently referred to with the IDs of the individual genes, joined by a plus (‘+’) sign. If you do not like this behaviour, you can disable aggregation with the option `-r no`. Without aggregation, exons that overlap with other exons from different genes are simply skipped.

    python /home/parashar/anaconda3/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py -r no /home/parashar/archive/rnaseq/gencode.v33.annotation.gtf gencode.v33.gff 2> gff.stderr

DEXSeq doesnot work with BAM files though theyclaim it cam. So convert the bams to sam files using samtools view subcommand.

    while read p
    do
    name=$(basename ${p} .sorted.sam)
    bsub -o $name.o -e $name.e -n 10 -q regularq "python /home/parashar/anaconda3/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -s reverse -a 10 gencode.v33.gff $p ${name}.Readcounts.txt 2> ${name}.stderr"
    done < list

