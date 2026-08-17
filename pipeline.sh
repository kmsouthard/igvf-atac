#!/usr/bin/env bash
#
# Bulk ATAC-seq: FASTQ to filtered BAM and CPM-normalized bigWig.
#
#     ./pipeline.sh pipeline_info/samplesheet.valid.csv results/
#
set -euo pipefail

SAMPLESHEET=${1:-pipeline_info/samplesheet.valid.csv}
OUTDIR=${2:-results}

nextflow run nf-core/atacseq \
    -r 2.1.2 \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --genome GRCh38 \
    --read_length 100 \
    -profile singularity

#   $OUTDIR/bwa/merged_library/<SAMPLE>.mLb.clN.sorted.bam
#   $OUTDIR/bwa/merged_library/<SAMPLE>.mLb.clN.sorted.bam.bai
#   $OUTDIR/bwa/merged_library/bigwig/<SAMPLE>.mLb.clN.bigWig
#   $OUTDIR/bwa/merged_library/bigwig/scale/<SAMPLE>.mLb.clN.scale_factor.txt
