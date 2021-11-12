# Generate synthetic reads
./hisat2_simulate_reads.py 22.fa 22.gtf 22.snp test -n 10000

# Build hisat2-quant binary program
make hisat2-quant-bin

# Run hisat-quant-bin
./hisat2-quant-bin


# Directories
#   hisat2/evaluation/tests/allele_specific

