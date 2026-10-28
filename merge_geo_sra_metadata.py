
#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

geo_path = Path(
    "data/processed_metadata/GSE178994_sample_metadata.csv"
)
sra_path = Path(
    "data/raw_metadata/sra/GSE178994_sra_metadata.tsv"
)
output_path = Path(
    "data/processed_metadata/GSE178994_merged_metadata.csv"
)

with geo_path.open(newline="", encoding="utf-8") as f:
    geo_reader = csv.DictReader(f)
    geo_fields = geo_reader.fieldnames
    geo_rows = list(geo_reader)

with sra_path.open(newline="", encoding="utf-8") as f:
    sra_rows = list(csv.DictReader(f, delimiter="\t"))

by_srx = defaultdict(list)

for row in sra_rows:
    if row.get("experiment_accession"):
        by_srx[row["experiment_accession"]].append(row)

sra_mapping = {
    "SRS": "sample_accession",
    "SRR": "run_accession",
    "SRP": "study_accession",
    "BioProject_SRA": "bioproject",
    "library_layout_SRA": "library_layout",
    "run_total_spots": "run_total_spots",
    "run_total_bases": "run_total_bases",
    "instrument_model_SRA": "instrument_model",
}

merged = []
unmatched = []
mismatches = []

for sample in geo_rows:
    gsm = sample["GSM"]
    matches = by_srx.get(sample["SRX"], [])

    if not matches:
        unmatched.append(gsm)
        continue

    for run in matches:
        geo_biosample = sample.get("BioSample", "")
        sra_biosample = run.get("biosample", "")

        if (
            geo_biosample
            and sra_biosample
            and geo_biosample != sra_biosample
        ):
            mismatches.append(gsm)
            continue

        record = dict(sample)

        for output_field, sra_field in sra_mapping.items():
            record[output_field] = run.get(sra_field, "")

        merged.append(record)

if unmatched or mismatches:
    raise ValueError(
        f"Unmatched GSMs: {unmatched}; "
        f"BioSample mismatches: {mismatches}"
    )

fields = list(geo_fields)

for field in sra_mapping:
    if field not in fields:
        fields.append(field)

output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(merged)

print("GEO samples:", len(geo_rows))
print("SRA runs:", len(sra_rows))
print("Merged rows:", len(merged))
print("Unmatched samples:", len(unmatched))
print("BioSample mismatches:", len(mismatches))
print("Output:", output_path)