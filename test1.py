import csv
from collections import Counter

path = "data/processed_metadata/GSE178994_sample_metadata.csv"

with open(path, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

groups = Counter(
    (r["disease_state"], r["stimulation"])
    for r in rows
)

print("Disease × Stimulation:")
for (disease, stimulation), count in sorted(groups.items()):
    print(f"{disease:20} | {stimulation:30} | {count}")

print("\nNon-terminal-ileum samples:")
for r in rows:
    tissue = r["tissue_source"].strip().casefold()
    if tissue != "terminal ileum":
        print(r["GSM"], r["tissue_source"],
              r["disease_state"], r["stimulation"])