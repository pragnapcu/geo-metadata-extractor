import csv
from collections import defaultdict

geo_file = "data/processed_metadata/GSE178994_sample_metadata.csv"
sra_file = "data/raw_metadata/sra/GSE178994_sra_metadata.tsv"

with open(geo_file, newline="", encoding="utf-8") as f:
    geo = list(csv.DictReader(f))

with open(sra_file, newline="", encoding="utf-8") as f:
    sra = list(csv.DictReader(f, delimiter="\t"))

by_srx = defaultdict(list)
for row in sra:
    by_srx[row["experiment_accession"]].append(row)

matched = 0
unmatched = []
biosample_mismatches = []

for sample in geo:
    gsm = sample["GSM"]
    srx = sample["SRX"]
    hits = by_srx.get(srx, [])

    if not hits:
        unmatched.append(gsm)
        continue

    matched += 1

    for run in hits:
        if sample["BioSample"] != run["biosample"]:
            biosample_mismatches.append(gsm)

        print(
            gsm,
            run["sample_accession"],
            srx,
            run["run_accession"],
            sep="\t"
        )

print("\nGEO samples:", len(geo))
print("SRA runs:", len(sra))
print("Matched GEO samples:", matched)
print("Unmatched GEO samples:", unmatched)
print("BioSample mismatches:", biosample_mismatches)