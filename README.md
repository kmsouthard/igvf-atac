# igvf-atac

Bulk ATAC-seq processing for Hs27 fibroblasts and hTERT RPE-1 cells, from
*Comprehensive transcription factor perturbations recapitulate fibroblast
transcriptional states* (Nat Genet 2025,
[doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)).

## Run

```bash
./pipeline.sh pipeline_info/samplesheet.valid.csv results/
```

Requires Nextflow 23.10.0 and Singularity. All other software is pulled by the
pipeline; `pipeline_info/software_versions.yml` records the versions of the run that
produced the deposited files.

## Inputs

`pipeline_info/samplesheet.valid.csv` — 6 Hs27 and 3 RPE-1 libraries. Paths in it are
the raw FASTQs, available from SRA
[PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254).

`--genome GRCh38` resolves through iGenomes to the reference FASTA, BWA index and
ENCODE blacklist. Nothing else is needed locally.

## Outputs

Per replicate, under `results/bwa/merged_library/`:

```
<SAMPLE>.mLb.clN.sorted.bam        merged library, filtered alignments
<SAMPLE>.mLb.clN.sorted.bam.bai
bigwig/<SAMPLE>.mLb.clN.bigWig     coverage, counts per million
bigwig/scale/<SAMPLE>.mLb.clN.scale_factor.txt
```

`.clN` is MAPQ-filtered, duplicate-removed, chrM-dropped, blacklist-filtered and
orphan-stripped. The bigWig scale factor is `1e6 / mapped reads`.

Peak calling, consensus peaks, featureCounts, DESeq2, ataqv and the pooled `.mRp.clN`
merged-replicate outputs are also produced under `results/`.

## Data

IGVF Data Portal, analysis set
[IGVFDS5661FPYR](https://data.igvf.org/analysis-sets/IGVFDS5661FPYR/) — the three Hs27
`.mLb.clN` BAMs, indexes and bigWigs. Its inputs are measurement sets
[IGVFDS4038GSGB](https://data.igvf.org/measurement-sets/IGVFDS4038GSGB/),
[IGVFDS7678SSQX](https://data.igvf.org/measurement-sets/IGVFDS7678SSQX/) and
[IGVFDS8072MCFE](https://data.igvf.org/measurement-sets/IGVFDS8072MCFE/).

These records are `in progress`, so the links resolve only for signed-in submitters
until release.

## Related

* Pipeline: [nf-core/atacseq](https://nf-co.re/atacseq/2.1.2) 2.1.2
* [kmsouthard/igvf-cutrun](https://github.com/kmsouthard/igvf-cutrun) ·
  [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra)
