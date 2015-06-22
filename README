# ChIP-seq Analysis Usage:

# 0) Initialization
####################################################################
SPECIES=crypto
GENE=USV1
CONDITION=inducing
LIBRARY=H99/bowtie_index/crNeoH99
GENOMIC_FEATURE=H99/crNeoH99.gff3

# 1) Align the reads to the genome
####################################################################
python build_stages.py -stage 1 -species $SPECIES -gene $GENE -condition $CONDITION -bowtie_library $LIBRARY
job_name=$SPECIES.$GENE.$CONDITION.rBowtie
nq $SPECIES/$GENE/$CONDITION/analysis/rBowtie | qsub -P long -l h_vmem=4G -N ${job_name}

# 2) Identify significant peaks of aligned reads
#################################################################### 
python build_stages.py -stage 2 -species $SPECIES -gene $GENE -condition $CONDITION
job_name=$SPECIES.$GENE.$CONDITION.rMACS
for x in $SPECIES/$GENE/$CONDITION/analysis/peaks_*; do qsub -V -cwd -N ${job_name} $x/rMACS; done

# 3) Identify putative bound genes
####################################################################
python build_stages.py -stage 3 -species $SPECIES -gene $GENE -condition $CONDITION -genomic_feature=${GENOMIC_FEATURE}
for x in $SPECIES/$GENE/$CONDITION/analysis/peaks_*; do bash $x/rAnnotatePeaks; done

# 4) Clean up log files
####################################################################
mkdir -p log
mv $SPECIES.$GENE.$CONDITION.* log
