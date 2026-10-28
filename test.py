import csv
from collections import Counter

path = "data/processed_metadata/GSE178994_sample_metadata.csv"

with open(path, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

print("Total samples:", len(rows))

for field in ["disease_state", "stimulation", "tissue_source"]:
    print(f"\n{field}:")
    for value, count in Counter(row.get(field, "") for row in rows).items():
        print(f"  {value}: {count}")

print("\nFirst sample:")
for field in ["GSM", "sample_title", "tissue_source",
              "cell_type", "disease_state", "stimulation",
              "BioSample", "SRX"]:
    print(f"{field}: {rows[0].get(field, '')}")