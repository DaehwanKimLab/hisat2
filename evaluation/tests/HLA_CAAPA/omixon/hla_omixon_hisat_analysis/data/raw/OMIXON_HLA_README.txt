HLA raw data:

For each locus, the "Winner" calls are based on %of locus covered and average read coverage. The first pair (or allele if the second allele is Identical) is the "Winner" and is the call for that locus for that individual. Other close calls are also given, but winner: false for these calls.

For some loci, these coverage stats are used to classify the calls into good or low quality based on where allele-defining variants are known to occur.

These quality control measures are noted in the HLA output table of the types for all individuals for the loci: A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1.The numbers in brackets were used for QC checks. The first two numbers correspond to exon and the second two numbers correspond to exon 3. The first number in each pair is % coverage, or how much of the exon is covered by the reads. The second number in each pair is Average Overlap, basically the average read coverage for the exon. Depending on which sites are known to describe different alleles for a locus, either both exon2 and exon3 had to pass QC (% coverage >= 85 and average read coverage >= 5) or only exon2 had to pass QC. For A,B,C, DRB1, DQA1 exons 2 & 3 were used. For DQB1, DPA1, DPB1 only exon 2 was used. We worked with Omixon to generate these criteria. They said no such criteria have been established for other loci.

