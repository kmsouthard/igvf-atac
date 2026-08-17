#!/usr/bin/env python3
"""Check the Hs27 bulk ATAC files and the replicate mapping against the nf-core run.

Two jobs. The first is the usual one: confirm the submitted files are the run's
output and that the normalization claim holds. The second is the one that matters
here -- nf-core numbered its samples independently of the bench, so this checks
which library is actually inside each output file, three separate ways, and
compares that against what the portal says each file derives from.

Dependency-free. Exits non-zero if any check fails.

    python3 verify/verify.py <nf-core results dir> [--md5]

--md5 also hashes the BAMs, which are 2-4 GB each and slow.
"""
import csv
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bw_quantum as bq
from bw_summary import read_bw

CHR1_GRCH38 = 248_956_422
CPM_TARGET = 1_000_000

# accession, md5 and the declared derived_from are from the IGVF portal.
# library and fastq_pairs are what the nf-core samplesheet says is inside.
FILES = [
    dict(nfcore="Hs27_REP1", library="RA479-3", pairs=1,
         bam=("IGVFFI8734PFCK", "1bd4288927f1790e8cbed88c6fdb6509"),
         bai=("IGVFFI1130AVUX", "a3a10ca70030a84f308f10dfe71b8790"),
         bigwig=("IGVFFI4199OOUU", "30c42697a7517795c5e18f71565e8e4c"),
         portal_derived_from=4),
    dict(nfcore="Hs27_REP2", library="RA479-2", pairs=3,
         bam=("IGVFFI1333TQNF", "62c38486d05b38fc73ed78928fb905ec"),
         bai=("IGVFFI4980HZQV", "1ba737907be9632f00429d7b5d56aa91"),
         bigwig=("IGVFFI0367LVLH", "10b64218fe790d920c42668cbe56834c"),
         portal_derived_from=6),
    dict(nfcore="Hs27_REP3", library="RA479-1", pairs=2,
         bam=("IGVFFI8359FWPE", "c09c8098e05407b78633e010d5b09f8d"),
         bai=("IGVFFI6565PXAY", "a14a02dbf0c5de633b3d37990ac8de81"),
         bigwig=("IGVFFI8794JKNW", "c30b6918fa5e3bf76a5f361cabe0c922"),
         portal_derived_from=2),
]


class Report:
    def __init__(self):
        self.passed = self.failed = 0

    def check(self, cond, what, detail=""):
        if cond:
            self.passed += 1
            print(f"  PASS  {what}" + (f"  [{detail}]" if detail else ""))
        else:
            self.failed += 1
            print(f"  FAIL  {what}" + (f"  [{detail}]" if detail else ""))
        return cond


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def samplesheet_libraries(results):
    """nf-core sample name -> (library, number of fastq pairs), from the samplesheet."""
    path = os.path.join(results, "pipeline_info", "samplesheet.valid.csv")
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            sample = row["sample"].rsplit("_T", 1)[0]
            lib = os.path.basename(row["fastq_1"]).split("_IGO")[0]
            libs, n = out.get(sample, (set(), 0))
            libs.add(lib)
            out[sample] = (libs, n + 1)
    return {k: (sorted(v[0]), v[1]) for k, v in out.items()}


def mapped_reads(results, sample):
    path = os.path.join(results, "bwa", "merged_library", "samtools_stats",
                        f"{sample}.mLb.clN.sorted.bam.flagstat")
    with open(path) as fh:
        return int(fh.readline().split()[0])


def scale_factor(results, sample):
    path = os.path.join(results, "bwa", "merged_library", "bigwig", "scale",
                        f"{sample}.mLb.clN.scale_factor.txt")
    with open(path) as fh:
        return float(fh.read().strip())


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    results = sys.argv[1]
    do_bam = "--md5" in sys.argv

    sheet = samplesheet_libraries(results)
    r = Report()
    ml = os.path.join(results, "bwa", "merged_library")

    print("\n=== which library is in which output file ===\n")
    for f in FILES:
        s = f["nfcore"]
        libs, pairs = sheet.get(s, ([], 0))
        print(f"{s}  portal says derived from {f['portal_derived_from']} fastqs")

        r.check(libs == [f["library"]],
                f"samplesheet: {s} was built from one library",
                f"{libs} x {pairs} pair(s)")
        r.check(pairs == f["pairs"], "pair count matches the crosswalk", f"{pairs}")
        r.check(f["portal_derived_from"] == 2 * pairs,
                "portal derived_from count matches the samplesheet",
                f"portal {f['portal_derived_from']}, nf-core {2 * pairs}")

        mapped = mapped_reads(results, s)
        sc = scale_factor(results, s)
        r.check(abs(CPM_TARGET / sc - mapped) / mapped < 1e-4,
                "scale factor is counts per million of mapped reads",
                f"{sc} vs 1e6/{mapped:,} = {CPM_TARGET / mapped:.7f}")

        bw = os.path.join(ml, "bigwig", f"{s}.mLb.clN.bigWig")
        r.check(md5(bw) == f["bigwig"][1],
                f"bigWig is byte-identical to {f['bigwig'][0]} on the portal")

        d = read_bw(bw)
        sizes = dict(d["chroms"])
        r.check(sizes.get("chr1") == CHR1_GRCH38, "assembly is GRCh38",
                f"chr1 = {sizes.get('chr1'):,}")

        fh, endian, meta = bq._open(bw)
        with fh:
            items = bq.leaves(fh, endian, meta["index"])
            step = max(1, len(items) // 200)
            vals = []
            for off, size in items[::step][:200]:
                try:
                    vals.extend(bq.block_values(fh, endian, meta["uncompress"], off, size))
                except Exception:                     # noqa: BLE001
                    continue
        q = bq.quantum(vals)
        r.check(abs(q - sc) / sc < 1e-4,
                "scale factor recovered from the bigWig matches the pipeline's own",
                f"recovered {q:.7f}, nf-core wrote {sc}")

        idx = os.path.join(ml, f"{s}.mLb.clN.sorted.bam.bai")
        if os.path.exists(idx):
            r.check(md5(idx) == f["bai"][1],
                    f"index is byte-identical to {f['bai'][0]} on the portal")
        if do_bam:
            bam = os.path.join(ml, f"{s}.mLb.clN.sorted.bam")
            r.check(md5(bam) == f["bam"][1],
                    f"BAM is byte-identical to {f['bam'][0]} on the portal")
        print()

    print("=" * 78)
    print(f"{'nf-core':10}{'contains':10}{'pairs':>7}{'mapped':>13}{'portal df':>11}  verdict")
    for f in FILES:
        libs, pairs = sheet.get(f["nfcore"], ([], 0))
        ok = f["portal_derived_from"] == 2 * pairs
        print(f"{f['nfcore']:10}{f['library']:10}{pairs:>7}"
              f"{mapped_reads(results, f['nfcore']):>13,}{f['portal_derived_from']:>11}"
              f"  {'ok' if ok else 'MIS-MAPPED'}")

    print("\n" + "=" * 78)
    print(f"{r.passed} passed, {r.failed} failed")
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
