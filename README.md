# igvf-atac

Processing pipeline for the Hs27 fibroblast bulk ATAC-seq dataset submitted to the
IGVF data portal, from *Comprehensive transcription factor perturbations recapitulate
fibroblast transcriptional states* (Nat Genet 2025,
[doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)).

Three biological replicates of parental Hs27 fibroblasts, paired-end 100 bp on a
NovaSeq 6000, processed with **nf-core/atacseq 2.1.2**. Nine submitted files:
a filtered BAM, its index, and a CPM-normalized bigWig per replicate, in analysis set
[IGVFDS5661FPYR](https://data.igvf.org/analysis-sets/IGVFDS5661FPYR/).

Unlike the [CUT&RUN companion](https://github.com/kmsouthard/igvf-cutrun), this run
survives intact — every submitted file is byte-identical to the output still on disk,
and `pipeline_info/` here is the run's own samplesheet and software manifest, not a
reconstruction.

## Read this before using the files

**The nf-core replicate numbers are not the bench replicate numbers.** nf-core assigned
`REP1`/`REP2`/`REP3` in samplesheet order; the bench numbered libraries `RA479-1/2/3`
independently. They do not correspond:

| File is named | Actually contains | Bench replicate | Measurement set | Fastq pairs | Mapped reads |
|---|---|---|---|---|---|
| `Hs27_REP1` | **`RA479-3`** | 3 | [IGVFDS8072MCFE](https://data.igvf.org/measurement-sets/IGVFDS8072MCFE/) | 1 | 43,527,858 |
| `Hs27_REP2` | `RA479-2` | 2 | [IGVFDS7678SSQX](https://data.igvf.org/measurement-sets/IGVFDS7678SSQX/) | 3 | 79,000,690 |
| `Hs27_REP3` | **`RA479-1`** | 1 | [IGVFDS4038GSGB](https://data.igvf.org/measurement-sets/IGVFDS4038GSGB/) | 2 | 74,051,896 |

Never link an analysis file by the replicate number in its filename. `results/replicate_crosswalk.csv`
is the mapping; `verify/verify.py` re-derives it three independent ways and checks it
against what the portal currently declares.

## The pipeline

`pipeline.sh` runs it; this is the command the submitted files came from, not a
reconstruction.

```bash
nextflow run nf-core/atacseq -r 2.1.2 \
    --input pipeline_info/samplesheet.valid.csv --outdir results \
    --genome GRCh38 --read_length 100 -profile singularity
```

nf-core/atacseq **2.1.2** on Nextflow 23.10.0, broad peak mode. The exact samplesheet
and the full software manifest are in `pipeline_info/`. Steps that produced the
submitted files:

| Step | Tool | Version |
|---|---|---|
| Read trimming | Trim Galore! / cutadapt / FastQC | 0.6.7 / 3.4 / 0.11.9 |
| Alignment | BWA · samtools | 0.7.17-r1188 · 1.16.1 |
| Filtering (`.clN`) | Picard MarkDuplicates · BamTools · samtools · BEDTools | 3.0.0 · 2.5.2 · 1.15.1 · 2.30.0 |
| Signal generation | BEDTools · UCSC bedGraphToBigWig | 2.30.0 · 445 |
| Sort and index | samtools | 1.17 |

**`.mLb.clN` is the merged-library cleaned alignment.** MAPQ-filtered, duplicates
removed by Picard, mitochondrial reads dropped, ENCODE blacklist regions removed, and
orphan mates stripped so only proper pairs survive. It is a heavily filtered BAM, which
is why `filtered: true` is the correct portal value for it.

**The bigWigs are counts per million**, not raw depth: `bedtools genomecov` is given a
per-file scale factor of `1e6 / mapped reads`, and nf-core writes those factors out
next to the tracks. Coverage is sparse — only contigs and positions with signal are
written, so contig counts differ per file (146, 158, 151) at 41–52% genome breadth.

Peak calling (MACS2 2.2.7.1), consensus peaks, featureCounts 2.0.1, DESeq2 1.28.0,
ataqv 1.3.1 and the pooled `mRp.clN` merged-replicate tracks all ran but are **not**
part of the IGVF submission. `metadata/analysis_step_version.csv` marks which steps are
needed and which are not.

## Verifying

```
python3 verify/verify.py /path/to/nextflow_atac/results        # bigWigs and indexes
python3 verify/verify.py /path/to/nextflow_atac/results --md5  # also hashes the BAMs
```

Per replicate it checks: which library the samplesheet says went in, that the fastq
pair count agrees, that the portal's `derived_from` count agrees with both, that the
scale factor is exactly `1e6 / mapped reads` from the flagstat, that the same factor
can be recovered from the bigWig's own value quantization, that the file is
byte-identical to its portal md5, and that the assembly is GRCh38.

**It currently exits non-zero: 22 passed, 2 failed.** Both failures are the
`derived_from` mis-mapping on `Hs27_REP1` and `Hs27_REP3` described above — a live
defect in the portal records, not in this repository. `metadata/patch_files.csv` has
the correction; the harness will pass once it is applied.

## Layout

```
pipeline.sh     the pipeline, end to end
pipeline_info/  the run's own samplesheet and software_versions.yml
verify/         verify.py, plus three standalone bigWig readers
results/        the replicate crosswalk
metadata/       the portal records these files imply, as submission sheets
```

## Related

* Paper: [doi:10.1038/s41588-025-02284-1](https://doi.org/10.1038/s41588-025-02284-1)
  · preprint [PMC11312553](https://pmc.ncbi.nlm.nih.gov/articles/PMC11312553/)
* Pipeline: [nf-core/atacseq](https://nf-co.re/atacseq/2.1.2) 2.1.2
* Raw reads: SRA [PRJNA1108254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1108254)
  · a separately-processed, 10M-normalized version of these libraries is on Zenodo
  [10.5281/zenodo.15215216](https://doi.org/10.5281/zenodo.15215216) and is *not* what
  was submitted here
* IGVF portal: [Tom Norman lab, MSKCC](https://data.igvf.org/labs/tom-norman/)
* Companion repositories:
  [kmsouthard/igvf-cutrun](https://github.com/kmsouthard/igvf-cutrun), `igvf-rnaseq`,
  and [kmsouthard/igvf-crispra](https://github.com/kmsouthard/igvf-crispra)
