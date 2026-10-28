import csv
from collections import Counter

path = "data/raw_metadata/sra/GSE178994_sra_metadata.tsv"

with open(path, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

print("SRA rows:", len(rows))
print("Columns:", list(rows[0]))

for field in ("run_accession", "experiment_accession",
              "sample_accession", "biosample"):
    counts = Counter(r.get(field, "") for r in rows)
    print(
        field,
        "unique:", len(counts),
        "duplicates:", sum(n - 1 for n in counts.values() if n > 1)
    )

print("\nGSM to SRR mapping:")
for r in rows:
    print(
        r.get("sample_accession", ""),
        r.get("experiment_accession", ""),
        r.get("run_accession", "")
    )