#!/usr/bin/env bash
#
# Hs27 bulk ATAC-seq: FASTQ to filtered BAM and CPM-normalized bigWig.
#
# This is the command that produced the submitted files, taken from the run's
# own execution report. Unlike the CUT&RUN companion this is not a
# reconstruction -- every submitted file is byte-identical to what this
# produced, which verify/verify.py checks.
#
#     ./pipeline.sh pipeline_info/samplesheet.valid.csv results/
#
set -euo pipefail

SAMPLESHEET=${1:-pipeline_info/samplesheet.valid.csv}
OUTDIR=${2:-results}

# nf-core/atacseq 2.1.2, Nextflow 23.10.0. --genome GRCh38 resolves through
# igenomes to a UCSC-named GRCh38 FASTA, its BWA index and the ENCODE blacklist.
# --read_length 100 selects the matching MACS2 effective genome size.
# Broad peak mode is the pipeline default and was not overridden.
nextflow run nf-core/atacseq \
    -r 2.1.2 \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    --genome GRCh38 \
    --read_length 100 \
    -profile singularity

# Outputs that were submitted to IGVF, per replicate:
#
#   $OUTDIR/bwa/merged_library/<SAMPLE>.mLb.clN.sorted.bam
#   $OUTDIR/bwa/merged_library/<SAMPLE>.mLb.clN.sorted.bam.bai
#   $OUTDIR/bwa/merged_library/bigwig/<SAMPLE>.mLb.clN.bigWig
#
# .mLb  merged library -- the per-replicate merge of its sequencing runs
# .clN  cleaned -- MAPQ filtered, duplicates removed, chrM dropped, ENCODE
#       blacklist regions removed, orphan mates stripped
#
# The bigWig is bedtools genomecov given --scale 1e6/mapped_reads, so it is
# counts per million rather than raw depth. The factor used for each file is
# written to bigwig/scale/<SAMPLE>.mLb.clN.scale_factor.txt and can also be
# recovered from the bigWig itself with verify/bw_quantum.py.
#
# Everything else the run produced -- MACS2 broad peaks, consensus peaks,
# featureCounts, DESeq2 QC, ataqv, and the pooled .mRp.clN merged-replicate
# tracks -- exists but was not submitted. See metadata/analysis_step_version.csv.
#
# NOTE ON SAMPLE NAMES: nf-core numbers replicates in samplesheet order, which
# is not the bench numbering. Hs27_REP1 contains RA479-3 and Hs27_REP3 contains
# RA479-1. See results/replicate_crosswalk.csv before linking anything.
