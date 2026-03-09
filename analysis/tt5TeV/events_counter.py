import os
import gzip
import pickle
import argparse
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# -----------------------------
# configuration
# -----------------------------
process_map = {
    "schan": ["schannel","sbarchannel"],
    "tt"   : ["ttPS"],
    "tW"   : ["tW","tbarW"],
    "tchan": ["tchannel","tbarchannel"],
    "DY"   : ["DYJetsToLLM10to50","DYJetsToLLMLL50"],
    "Wjets": ["W0JetsToLNu","W1JetsToLNu","W2JetsToLNu","W3JetsToLNu"],
    "QCD"  : ["QCD"],
    "data" : ["HighEGJet","SingleMuon"]
}

process_order = ["schan","tt","tW","tchan","DY","Wjets","QCD"]

histname = "counts"
syst = "norm"
lumi = 302

# -----------------------------
# compute yields
# -----------------------------
def compute_yields(folder, channel, level):
    counts = defaultdict(float)
    for file in os.listdir(folder):
        if not file.endswith(".pkl.gz"):
            continue
        path = os.path.join(folder, file)
        with gzip.open(path,"rb") as f:
            h = pickle.load(f)
        hist = h[histname]
        values = hist.sum(histname).values()
        for key,val in values.items():
            dataset, ch, lev, sy = key
            if ch!=channel or lev!=level or sy!=syst:
                continue
            for proc, datasets in process_map.items():
                if dataset in datasets:
                    counts[proc] += val
    return counts

# -----------------------------
# build PDF table
# -----------------------------
# -----------------------------
# build PDF table (with separate data row)
# -----------------------------
def save_pdf_table(results, regions, output="yields.pdf"):
    columns = ["Process"] + regions

    rows = []
    # Add MC process rows (exclude 'data')
    for proc in process_order:
        if proc == "data":
            continue
        row = [proc]
        for r in regions:
            val = results[r][proc] * lumi
            row.append(val)
        rows.append(row)

    # Add MC total row
    mc_total_row = ["MC total"]
    for r in regions:
        total = sum(results[r][p]*lumi for p in process_order if p != "data")
        mc_total_row.append(total)
    rows.append(mc_total_row)

    # Add data row
    data_row = ["data"]
    for r in regions:
        data_val = results[r]["data"]
        data_row.append(data_val)
    rows.append(data_row)

    # Add Data/MC row
    datamc_row = ["Data/MC"]
    for idx,r in enumerate(regions):
        mc = mc_total_row[idx+1]  # skip first column
        data_val = data_row[idx+1]
        ratio = data_val/mc if mc>0 else 0
        datamc_row.append(ratio)
    rows.append(datamc_row)

    # Create PDF
    fig, ax = plt.subplots(figsize=(1.5*len(columns),1.2*len(rows)))
    ax.axis("off")
    table = ax.table(
        cellText=[[f"{x:.2f}" if isinstance(x,float) else str(x) for x in row] for row in rows],
        colLabels=columns,
        loc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1,2)

    pp = PdfPages(output)
    pp.savefig(fig, bbox_inches="tight")
    pp.close()
    plt.close(fig)
    print(f"PDF table saved to {output}")

# -----------------------------
# main
# -----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input", required=True, help="input folder with pkls")
    parser.add_argument("-c","--channels", nargs="+", default=["m"], help="channels")
    parser.add_argument("-l","--levels", nargs="+", default=["2j2b"], help="levels")
    parser.add_argument("-o","--output", default="yields.pdf", help="PDF output file")
    args = parser.parse_args()

    results = {}
    regions = []

    for ch in args.channels:
        for lev in args.levels:
            region = f"{ch}_{lev}"
            print("Processing", region)
            counts = compute_yields(args.input, ch, lev)
            results[region] = counts
            regions.append(region)

    save_pdf_table(results, regions, output=args.output)

if __name__ == "__main__":
    main()
