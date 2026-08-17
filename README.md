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

## Note on sample names

nf-core numbers replicates in samplesheet order, which is not the bench numbering:
`Hs27_REP1` is library `RA479-3`, `Hs27_REP2` is `RA479-2`, `Hs27_REP3` is `RA479-1`.

## Related

* Pipeline: [nf-core/atacseq](https://nf-co.re/atacseq/2.1.2) 2.1.2
* [kmsouthard/igvf-cutrun](https://github.com/kmsouthard/igvf-cutrun) ·
  [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra)
