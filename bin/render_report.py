#!/usr/bin/env python3
"""
bin/render_report.py

"""
import argparse, csv, gzip, glob, json, random, sys
from datetime import datetime
from pathlib import Path

random.seed(42)

def load_qc():
    rows = []
    for f in sorted(glob.glob("*_qc_metrics.tsv")):
        try:
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    rows.append(row)
        except Exception as e:
            print(f"  [WARN] qc {f}: {e}", file=sys.stderr)
    return rows

def load_summary():
    rows = []
    for f in sorted(glob.glob("*_summary.tsv")):
        try:
            with open(f) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    rows.append(row)
        except Exception as e:
            print(f"  [WARN] summary {f}: {e}", file=sys.stderr)
    return rows

def load_sumstats():
    hits, manhattan = [], []
    for gz in sorted(glob.glob("*_sumstats.tsv.gz")):
        try:
            with gzip.open(gz, "rt") as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    try:
                        p  = float(row.get("P",  "1"))
                        lp = float(row.get("LOG10P", "0"))
                        if p < 5e-8 and len(hits) < 5000:
                            hits.append({k: row.get(k,"") for k in
                                ["CHR","BP","SNP","P","LOG10P","BETA","SE",
                                 "A1","A2","EAF","N","TRAIT","COHORT"]})
                        # Manhattan: all significant + 4% random sample
                        if p < 5e-5 or random.random() < 0.04:
                            manhattan.append({
                                "c":  row.get("CHR",""),
                                "b":  row.get("BP",""),
                                "s":  row.get("SNP",""),
                                "lp": row.get("LOG10P",""),
                                "be": row.get("BETA",""),
                                "t":  row.get("TRAIT",""),
                                "co": row.get("COHORT",""),
                            })
                    except: pass
        except Exception as e:
            print(f"  [WARN] {gz}: {e}", file=sys.stderr)
    print(f"[REPORT] GW-sig: {len(hits)}  Manhattan pts: {len(manhattan)}")
    return hits, manhattan

def load_pca(eigenvec_path, population_file):
    pop_map = {}
    if population_file and Path(population_file).exists():
        try:
            with open(population_file) as f:
                reader = csv.DictReader(f, delimiter="\t")
                cols = reader.fieldnames or []
                id_col   = next((c for c in cols
                                 if c.strip("#").upper() in ("IID","SAMPLE","ID")), None)
                pop_col  = next((c for c in cols
                                 if "POP" in c.upper() and "SUPER" not in c.upper()), None)
                spop_col = next((c for c in cols
                                 if "SUPER" in c.upper() or c.upper()=="SUPERPOP"), None)
                if id_col and pop_col:
                    for row in reader:
                        pop_map[row[id_col]] = {
                            "pop":  row.get(pop_col, ""),
                            "spop": row.get(spop_col, "") if spop_col else "",
                        }
            print(f"[REPORT] Population map: {len(pop_map)} samples")
        except Exception as e:
            print(f"  [WARN] population file: {e}", file=sys.stderr)

    pca_rows = []
    try:
        lines = Path(eigenvec_path).read_text().splitlines()
        for line in lines[1:]:
            p = line.split()
            if len(p) < 4: continue
            iid = p[1]
            info = pop_map.get(iid, {})
            spop = info.get("spop", "")
            # Derive cohort from super-population (not IID prefix)
            # Cohort_A = EUR + EAS, Cohort_B = AFR + AMR + SAS
            if spop in ("EUR", "EAS"):
                coh = "Cohort_A"
            elif spop in ("AFR", "AMR", "SAS"):
                coh = "Cohort_B"
            else:
                coh = "Unknown"
            pca_rows.append({
                "iid": iid,
                "pc1": p[2], "pc2": p[3],
                "pc3": p[4] if len(p) > 4 else "0",
                "pc4": p[5] if len(p) > 5 else "0",
                "pop":  info.get("pop",  ""),
                "spop": spop,
                "coh":  coh,
            })
    except Exception as e:
        print(f"  [WARN] PCA eigenvec: {e}", file=sys.stderr)
    return pca_rows

def load_evals(pca_variance_path):
    rows = []
    try:
        with open(pca_variance_path) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    except Exception as e:
        print(f"  [WARN] pca_variance: {e}", file=sys.stderr)
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eigenvec")
    ap.add_argument("--pca_variance")
    ap.add_argument("--template",    required=True)
    ap.add_argument("--author",      default="Nadeem Khan")
    ap.add_argument("--affiliation", default="Bioinformatician, INRS-CAFSB")
    ap.add_argument("--institute",   default="INRS-Centre Armand-Frappier Sante-Biotechnologie, Laval, QC, Canada")
    ap.add_argument("--github",      default="github.com/nkhan119")
    ap.add_argument("--population_file", default="")
    ap.add_argument("--maf",      default="0.01")
    ap.add_argument("--hwe",      default="1e-6")
    ap.add_argument("--geno",     default="0.02")
    ap.add_argument("--mind",     default="0.02")
    ap.add_argument("--win",      default="1000")
    ap.add_argument("--r2",       default="0.1")
    ap.add_argument("--npcs",     default="20")
    ap.add_argument("--s1bsize",  default="1000")
    ap.add_argument("--s2bsize",  default="400")
    args = ap.parse_args()

    print(f"[REPORT] Author    : {args.author}")
    print(f"[REPORT] Institute : {args.institute}")

    qc_rows     = load_qc()
    gwas_summary= load_summary()
    hits, man   = load_sumstats()
    pca_rows    = load_pca(args.eigenvec, args.population_file) if args.eigenvec else []
    evals_rows  = load_evals(args.pca_variance) if args.pca_variance else []

    tmpl = Path(args.template).read_text()

    for k, v in {
        "__QC_DATA__"       : json.dumps(qc_rows),
        "__GWAS_SUMMARY__"  : json.dumps(gwas_summary),
        "__GWAS_HITS__"     : json.dumps(hits[:3000]),
        "__MANHATTAN__"     : json.dumps(man[:15000]),
        "__PCA_DATA__"      : json.dumps(pca_rows[:5000]),
        "__EVALS_DATA__"    : json.dumps(evals_rows),
        "__AUTHOR__"        : args.author,
        "__AFFILIATION__"   : args.affiliation,
        "__INSTITUTE__"     : args.institute,
        "__GITHUB__"        : args.github,
        "__DATE__"          : datetime.now().strftime("%B %d, %Y"),
        "__MAF__"           : args.maf,
        "__HWE__"           : args.hwe,
        "__GENO__"          : args.geno,
        "__MIND__"          : args.mind,
        "__WIN__"           : args.win,
        "__R2__"            : args.r2,
        "__NPCS__"          : args.npcs,
        "__S1BSIZE__"       : args.s1bsize,
        "__S2BSIZE__"       : args.s2bsize,
    }.items():
        tmpl = tmpl.replace(k, v)

    Path("GWAS_Report.html").write_text(tmpl)
    size_mb = Path("GWAS_Report.html").stat().st_size / 1e6
    print(f"[REPORT] ✓ GWAS_Report.html  ({size_mb:.1f} MB)")

if __name__ == "__main__":
    main()
