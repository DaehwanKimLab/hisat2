# Generate synthetic reads
./hisat2_simulate_reads.py 22.fa 22.gtf 22.snp test -n 10000

# Build hisat2-quant binary program
make hisat2-quant-bin

# Run hisat-quant (wrapper of hisat2-quant-bin)
./hisat2-quant evaluation/tests/allele_specific/sim.gene.sam

# Directories
#   hisat2/evaluation/tests/allele_specific

