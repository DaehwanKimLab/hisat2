# Clone HISAT2
> git clone https://github.com/DaehwanKimLab/hisat2
> cd hisat2
hisat2> git checkout dev-hisat-gt

# Build hisat2-quant binary program
hisat2> make hisat2-quant-bin


# Copy 22.fa, 22.gtf, and 22.snp
hisat2/evaluation/tests/allele_specific> scp nucleus:/project/bioinformatics/Kim_lab/shared/GTQ/22_data.tar.gz .


# Generate synthetic reads
hisat2/evaluation/tests/allele_specific> ../../../hisat2_simulate_reads.py 22.fa 22.gtf 22.snp test -n 10000 -a


# Run hisat-quant (wrapper of hisat2-quant-bin).
hisat2/evaluation/tests/allele_specific> ../../../hisat2-quant test.gene.sam


