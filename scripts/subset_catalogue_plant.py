#!/usr/bin/env python3
"""Cut one reactor's community out of the Campanaro biogas MAG catalogue.

The catalogue (PRJNA602310) is 1,401 MAGs from 134 digester metagenomes across
many plants. Scored whole it is a specificity benchmark; presented whole it is
not a digester. Campanaro et al. 2020 Additional file 8 gives the relative
abundance of every MAG in every experiment, so a single plant's community can
be reconstructed: the MAGs that are actually abundant in that reactor, with
their abundances.

Inputs
  Additional file 8   13068_2020_1679_MOESM8_ESM.xlsx from figshare (CC BY),
                      https://ndownloader.figshare.com/files/21991575
                      Sheet "Additional file 8": 1,635 MAGs x 91 experiments,
                      relative abundance in percent of mapped reads.
  catalogue dir       what scripts/fetch_biogas_catalogue.sh wrote:
                      proteomes/<accession>.faa, catalogue_index.tsv,
                      catalogue_taxonomy.tsv
  isolate map         accession -> MAG code, from the NCBI assembly 'isolate'
                      field. Fetched once from the NCBI datasets API and cached;
                      needs network, so run this on a login node the first time.

The join is exact: the paper's Bin Id ends in the MAG code (AS24abBPME_57) and
NCBI stores that code verbatim as the assembly's isolate name. 1,401 of the
1,635 codes are deposited; the rest are medium-quality bins that never went to
NCBI, and they are listed separately so the denominator stays visible.

Output (default <catalogue>/plants/<plant>/) mirrors the catalogue layout, so
    CATALOGUE=<out> SAMPLE="<plant>" sbatch scripts/score_catalogue.sbatch
scores it unchanged:
  proteomes/<accession>.faa   symlinks (or copies with --copy)
  catalogue_taxonomy.tsv      the matching rows, header comments kept
  manifest.tsv                accession, MAG code, paper taxonomy, abundance
  not_deposited.tsv           abundant MAGs with no genome at NCBI
The manifest is what aggregate_community.py --abundance reads.

Stdlib only: the xlsx is read as the zip of XML that it is.
"""
import argparse, datetime, json, os, re, shutil, sys, urllib.request, zipfile
import xml.etree.ElementTree as ET

FIGSHARE_URL = "https://ndownloader.figshare.com/files/21991575"
SHEET = "Additional file 8"
NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome"
BIOPROJECT = "PRJNA602310"
SUMMARY_COLS = {"AVERAGE", "MAX", "MIN", "STD"}   # trailing formula columns
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "p": "http://schemas.openxmlformats.org/package/2006/relationships"}


# ---------------------------------------------------------------- xlsx reader
def _col_index(ref):
    """'CP3' -> 93 (zero-based column)."""
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_sheet(xlsx_path, sheet_name):
    """Return the named worksheet as a list of rows (lists of str or None)."""
    with zipfile.ZipFile(xlsx_path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {r.get("Id"): r.get("Target") for r in rels}
        target = None
        for s in wb.find("m:sheets", NS):
            if s.get("name") == sheet_name:
                target = rid_to_target[s.get(f"{{{NS['r']}}}id")]
        if target is None:
            names = [s.get("name") for s in wb.find("m:sheets", NS)]
            sys.exit(f"!! sheet {sheet_name!r} not in {xlsx_path}; sheets: {names}")
        target = target if target.startswith("xl/") else "xl/" + target.lstrip("/")
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")):
                shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
        rows = []
        for row in ET.fromstring(z.read(target)).iter(f"{{{NS['m']}}}row"):
            cells = {}
            for c in row:
                v = c.find("m:v", NS)
                if v is None:
                    inl = c.find("m:is", NS)
                    val = "".join(t.text or "" for t in inl.iter(f"{{{NS['m']}}}t")) if inl is not None else None
                elif c.get("t") == "s":
                    val = shared[int(v.text)]
                else:
                    val = v.text          # numbers and cached formula results
                cells[_col_index(c.get("r"))] = val
            width = max(cells) + 1 if cells else 0
            rows.append([cells.get(i) for i in range(width)])
    return rows


def load_abundance(xlsx_path):
    """-> (experiments, [(bin_id, mag_code, mag_name, taxonomy, {exp: pct})])"""
    rows = read_sheet(xlsx_path, SHEET)
    header = None
    for i, r in enumerate(rows):
        if r and r[0] == "Bin Id":
            header = i
            break
    if header is None:
        sys.exit("!! no 'Bin Id' header row in the sheet")
    hdr = rows[header]
    exps, idx = [], []
    for j, name in enumerate(hdr[3:], start=3):
        if not name or name in SUMMARY_COLS or name.startswith("TIMES "):
            continue
        exps.append(name)
        idx.append(j)
    mags = []
    for r in rows[header + 1:]:
        if not r or not r[0]:
            continue
        m = re.search(r"(AS\d+\w+?_\d+)$", str(r[0]))
        code = m.group(1) if m else None
        ab = {}
        for name, j in zip(exps, idx):
            v = r[j] if j < len(r) else None
            try:
                ab[name] = float(v) if v not in (None, "") else 0.0
            except ValueError:
                ab[name] = 0.0
        mags.append((r[0], code, r[1] or "", r[2] or "", ab))
    return exps, mags


# --------------------------------------------------------------- NCBI isolates
def load_isolates(cache_path):
    """accession -> MAG code, from the cache or the NCBI datasets API."""
    if os.path.isfile(cache_path) and os.path.getsize(cache_path) > 0:
        out = {}
        for line in open(cache_path):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out[f[0]] = f[1]
        return out
    print(f"==> fetching isolate names for {BIOPROJECT} from NCBI (login node needed)",
          file=sys.stderr)
    out, token, page = {}, None, 0
    while True:
        url = f"{NCBI_API}/bioproject/{BIOPROJECT}/dataset_report?page_size=1000"
        if token:
            url += f"&page_token={token}"
        with urllib.request.urlopen(url, timeout=180) as fh:
            d = json.load(fh)
        for rep in d.get("reports", []):
            iso = ((rep.get("organism") or {}).get("infraspecific_names") or {}).get("isolate", "")
            if iso:
                out[rep["accession"]] = iso
        page += 1
        token = d.get("next_page_token")
        print(f"    page {page}: {len(out)} cumulative", file=sys.stderr)
        if not token:
            break
    with open(cache_path, "w") as fh:
        fh.write("#accession\tisolate\n")
        for a, i in sorted(out.items()):
            fh.write(f"{a}\t{i}\n")
    return out


def load_index(path):
    """catalogue_index.tsv -> set of the deduplicated accessions the fetch kept."""
    acc = set()
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f and f[0]:
            acc.add(f[0])
    return acc


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--xlsx", help="Additional file 8 (.xlsx). Default <catalogue>/campanaro_additional_file_8.xlsx")
    ap.add_argument("--download", action="store_true",
                    help=f"fetch the xlsx from figshare if missing ({FIGSHARE_URL})")
    ap.add_argument("--catalogue", default=os.path.join(os.environ.get("DB", ""),
                                                        "ad_panel_test", "biogas_catalogue"),
                    help="dir written by fetch_biogas_catalogue.sh (default $DB/ad_panel_test/biogas_catalogue)")
    ap.add_argument("--isolates", help="accession->isolate cache TSV (default <catalogue>/catalogue_isolates.tsv)")
    ap.add_argument("--list", action="store_true", help="list experiments with MAG counts and exit")
    ap.add_argument("--plant", help="experiment column to extract, e.g. BP-Gimenells")
    ap.add_argument("--min-abundance", type=float, default=0.1,
                    help="keep MAGs at or above this relative abundance in percent (default 0.1)")
    ap.add_argument("--out", help="output dir (default <catalogue>/plants/<plant>)")
    ap.add_argument("--copy", action="store_true", help="copy proteomes instead of symlinking")
    a = ap.parse_args()

    if not a.catalogue or not os.path.isdir(a.catalogue):
        sys.exit(f"!! catalogue dir not found: {a.catalogue!r} (set DB or --catalogue)")
    xlsx = a.xlsx or os.path.join(a.catalogue, "campanaro_additional_file_8.xlsx")
    if not os.path.isfile(xlsx):
        if not a.download:
            sys.exit(f"!! {xlsx} missing; pass --download or --xlsx")
        print(f"==> downloading Additional file 8 -> {xlsx}", file=sys.stderr)
        urllib.request.urlretrieve(FIGSHARE_URL, xlsx)

    exps, mags = load_abundance(xlsx)
    n_coded = sum(1 for m in mags if m[1])
    print(f"==> {len(mags)} MAGs x {len(exps)} experiments in {os.path.basename(xlsx)}; "
          f"{n_coded} Bin Ids carry a MAG code", file=sys.stderr)

    if a.list or not a.plant:
        print(f"{'experiment':<40} {'MAGs>=' + str(a.min_abundance) + '%':>12} {'>=1%':>6} {'sum%':>7}")
        for e in exps:
            vals = [m[4][e] for m in mags]
            n = sum(1 for v in vals if v >= a.min_abundance)
            n1 = sum(1 for v in vals if v >= 1.0)
            print(f"{e:<40} {n:>12} {n1:>6} {sum(vals):>7.1f}")
        if not a.plant:
            print("\nBP-* are full-scale plants, LSBR-* lab-scale reactors. Pick one with --plant.",
                  file=sys.stderr)
        return 0

    if a.plant not in exps:
        close = [e for e in exps if a.plant.lower() in e.lower()]
        sys.exit(f"!! {a.plant!r} is not an experiment column" +
                 (f"; did you mean {close}?" if close else "; use --list"))

    index_path = os.path.join(a.catalogue, "catalogue_index.tsv")
    tax_path = os.path.join(a.catalogue, "catalogue_taxonomy.tsv")
    prot_dir = os.path.join(a.catalogue, "proteomes")
    for p in (index_path, tax_path, prot_dir):
        if not os.path.exists(p):
            sys.exit(f"!! {p} missing; run scripts/fetch_biogas_catalogue.sh first")
    kept = load_index(index_path)
    isolates = load_isolates(a.isolates or os.path.join(a.catalogue, "catalogue_isolates.tsv"))
    # The fetch kept one accession per genome (RefSeq over GenBank); honour that
    # choice so a MAG deposited under both never appears twice.
    code_to_acc = {}
    for acc, code in isolates.items():
        if acc in kept:
            code_to_acc[code] = acc
    print(f"==> {len(isolates)} isolate names from NCBI, {len(code_to_acc)} map onto the "
          f"{len(kept)} accessions the fetch kept", file=sys.stderr)

    selected = [m for m in mags if m[4][a.plant] >= a.min_abundance]
    selected.sort(key=lambda m: -m[4][a.plant])
    total_all = sum(m[4][a.plant] for m in mags)
    total_sel = sum(m[4][a.plant] for m in selected)
    deposited = [m for m in selected if m[1] in code_to_acc]
    missing = [m for m in selected if m[1] not in code_to_acc]
    total_dep = sum(m[4][a.plant] for m in deposited)
    no_faa = [m for m in deposited
              if not os.path.isfile(os.path.join(prot_dir, code_to_acc[m[1]] + ".faa"))]
    if no_faa:
        print(f"!! {len(no_faa)} selected genomes have no proteome under {prot_dir}; "
              f"the fetch is incomplete: {[code_to_acc[m[1]] for m in no_faa][:5]}...",
              file=sys.stderr)
        deposited = [m for m in deposited if m not in no_faa]

    out = a.out or os.path.join(a.catalogue, "plants", a.plant)
    os.makedirs(os.path.join(out, "proteomes"), exist_ok=True)
    for m in deposited:
        acc = code_to_acc[m[1]]
        src = os.path.join(prot_dir, acc + ".faa")
        dst = os.path.join(out, "proteomes", acc + ".faa")
        if os.path.lexists(dst):
            os.remove(dst)
        if a.copy:
            shutil.copyfile(src, dst)
        else:
            os.symlink(os.path.relpath(src, os.path.dirname(dst)), dst)

    # taxonomy: the catalogue rows for these accessions, warning header intact
    want = {code_to_acc[m[1]] for m in deposited}
    n_tax = 0
    with open(tax_path) as fi, open(os.path.join(out, "catalogue_taxonomy.tsv"), "w") as fo:
        for line in fi:
            if line.startswith("#") or line.startswith("user_genome"):
                fo.write(line)
                continue
            if line.split("\t", 1)[0] in want:
                fo.write(line)
                n_tax += 1

    stamp = datetime.date.today().isoformat()
    prov = [
        f"# plant: {a.plant}",
        f"# source: Campanaro et al. 2020 Additional file 8 ({FIGSHARE_URL}), sheet '{SHEET}'",
        f"# cutoff: relative abundance >= {a.min_abundance} % of reads mapped to the 1,635-MAG catalogue",
        f"# extracted: {stamp}",
        f"# mags_in_experiment_total_pct: {total_all:.2f}",
        f"# mags_above_cutoff: {len(selected)}",
        f"# mags_above_cutoff_pct: {total_sel:.2f}",
        f"# mags_scored: {len(deposited)}",
        f"# mags_scored_pct: {total_dep:.2f}",
        f"# mags_not_deposited: {len(missing)}",
        f"# mags_not_deposited_pct: {total_sel - total_dep:.2f}",
        "# mags_not_deposited are listed in not_deposited.tsv",
        "# abundance is a share of mapped reads across all catalogue MAGs, not of the raw sample",
    ]
    with open(os.path.join(out, "manifest.tsv"), "w") as fh:
        fh.write("\n".join(prov) + "\n")
        fh.write("accession\tmag_code\tmag_name\tpaper_taxonomy\trel_abundance_pct\tbin_id\n")
        for m in deposited:
            fh.write(f"{code_to_acc[m[1]]}\t{m[1]}\t{m[2]}\t{m[3]}\t{m[4][a.plant]:.6g}\t{m[0]}\n")
    with open(os.path.join(out, "not_deposited.tsv"), "w") as fh:
        fh.write("# MAGs above the cutoff in this experiment with no genome under "
                 f"{BIOPROJECT} (medium-quality bins). Their abundance is part of the "
                 "community but cannot be scored.\n")
        fh.write("bin_id\tmag_code\tmag_name\tpaper_taxonomy\trel_abundance_pct\n")
        for m in missing:
            fh.write(f"{m[0]}\t{m[1] or ''}\t{m[2]}\t{m[3]}\t{m[4][a.plant]:.6g}\n")

    print(f"==> {a.plant}: {len(selected)} MAGs >= {a.min_abundance} % "
          f"({total_sel:.1f} % of mapped reads); {len(deposited)} deposited and linked "
          f"({total_dep:.1f} %); {len(missing)} not deposited; {n_tax} taxonomy rows",
          file=sys.stderr)
    meth = sum(1 for m in deposited if "Euryarchaeota" in m[3] or "Methano" in m[3]
               or "Halobacterota" in m[3] or "Thermoplasmatota" in m[3])
    print(f"    {meth} of the linked MAGs carry an archaeal/methanogen taxonomy label "
          f"(paper's assignment, not the panel's call)", file=sys.stderr)
    print(f"==> wrote {out}", file=sys.stderr)
    print(f"Next: CATALOGUE={out} SAMPLE=\"{a.plant}\" sbatch -p <partition> scripts/score_catalogue.sbatch",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
