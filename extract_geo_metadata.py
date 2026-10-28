#!/usr/bin/env python3
"""Extract sample metadata from a GEO family SOFT file; Python >= 3.8, stdlib only.

No network requests or expression-table processing are performed. Standardized
biological fields use channel 1 only, with explicit attribute-name aliases;
no values are inferred from titles, protocols, diagnoses, or other channels.

CSV conventions:
* Empty cell: absent standardized/dynamic value. A single value is plain text;
  repeated values are a JSON array, including identical repeated occurrences.
* characteristic_chN__*: dynamically discovered characteristic attributes.
* geo__*: every original Sample metadata field, always a JSON array of strings.
* geo_metadata_json: ordered [original field name, original value] pairs. This
  is the authoritative, lossless representation of parsed sample metadata.
* standardized_provenance_json: zero-based indices into geo_metadata_json.
  GSE associations instead reference series membership/sample series_id/sole
  family series, as recorded in gse_association_source.
* geo_unparsed_lines_json: [line number, text] for malformed sample lines.

Only line terminators and ONE optional formatting space after '=' are removed
from original values. Characteristic parsing splits on the FIRST colon only.
Normalized fields trim surrounding whitespace; original values remain intact.
Distinct labels that normalize to the same column name receive hash suffixes.
Data tables and their column descriptions are intentionally excluded.

Fatal errors (including duplicate GSM sections and count mismatches) leave an
existing output untouched. Other validation findings warn, or fail with --strict.
Exit codes: 0 success, 1 input/validation/output failure, 2 invalid CLI arguments.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys
import tempfile
import zlib
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple


STANDARD_FIELDS = (
    "GSE", "GSM", "sample_title", "tissue_source", "cell_type",
    "disease_state", "stimulation", "treatment", "organism", "extracted_molecule",
    "library_strategy", "library_source", "library_selection",
    "instrument_model", "platform_id", "treatment_protocol", "growth_protocol",
    "extraction_protocol", "data_processing", "BioSample", "SRX",
    "supplementary_files", "sample_description", "SRR", "SRP", "BioProject",
)
IMPORTANT_FIELDS = ("GSE", "sample_title", "tissue_source", "organism",
                    "library_strategy")
# A validation expectation supplied by this project's specification, not a
# biological inference. Other studies have no implicit expected sample count.
KNOWN_COUNTS = {"GSE178994": 14}
DIRECT_FIELDS = {
    "sample_title": "title", "organism": "organism_ch1",
    "extracted_molecule": "molecule_ch1",
    "library_strategy": "library_strategy", "library_source": "library_source",
    "library_selection": "library_selection", "instrument_model": "instrument_model",
    "platform_id": "platform_id", "sample_description": "description",
    "treatment_protocol": "treatment_protocol_ch1",
    "growth_protocol": "growth_protocol_ch1",
    "extraction_protocol": "extract_protocol_ch1", "data_processing": "data_processing",
}
ALIASES = {
    "tissue_source": {"tissue", "tissue_source", "tissue_type", "source_name"},
    "cell_type": {"cell_type"},
    "disease_state": {"disease", "disease_state", "disease_status", "diagnosis"},
    "stimulation": {"stimulation", "stimulus"},
    "treatment": {"treatment", "treatment_condition", "drug_treatment", "drug",
                  "medication", "therapy"},
}
ENTITY = re.compile(r"^\^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.*?)\s*$")
CHARACTERISTIC = re.compile(r"^characteristics_ch([1-9][0-9]*)$", re.I)
TABLE = re.compile(r"^!(?:[A-Za-z0-9_]+)_table_(begin|end)\s*$", re.I)
ACCESSIONS = {"BioSample": re.compile(r"\b(?:SAMN\d+|SAMEA\d+|SAMD\d+)\b"),
              "SRX": re.compile(r"\bSRX\d+\b"),
              "SRR": re.compile(r"\bSRR\d+\b"),
              "SRP": re.compile(r"\bSRP\d+\b"),
              "BioProject": re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b")}


class MetadataError(Exception):
    """A readable extraction/validation failure."""


@dataclass
class Sample:
    gsm: str
    line: int
    metadata: List[Tuple[str, str]] = field(default_factory=list)
    unparsed: List[Tuple[int, str]] = field(default_factory=list)


@dataclass
class Study:
    samples: List[Sample] = field(default_factory=list)
    series: List[str] = field(default_factory=list)
    membership: Dict[str, List[str]] = field(default_factory=lambda: defaultdict(list))
    warnings: List[str] = field(default_factory=list)


def slug(text: str) -> str:
    """Normalize labels, never metadata values."""
    return re.sub(r"[^a-z0-9]+", "_", text.casefold()).strip("_") or "unnamed"


def to_json(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def cell(values: Sequence[str]) -> str:
    """Keep multiplicity; original untrimmed values are also stored separately."""
    if not values:
        return ""
    return values[0] if len(values) == 1 else to_json(list(values))


def open_soft(path: Path, encoding: str) -> TextIO:
    """Detect gzip by magic bytes, also rejecting damaged .gz files clearly."""
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if path.suffix.lower() == ".gz" and not compressed:
        raise MetadataError("Input has .gz suffix but no gzip header")
    opener = gzip.open if compressed else open
    return opener(path, "rt", encoding=encoding, errors="strict", newline="")


def parse_soft(path: Path, encoding: str = "utf-8-sig") -> Study:
    """Stream SOFT entities; retain only sample metadata and series membership."""
    study = Study()
    current: Optional[Sample] = None
    entity_kind, series_id = "", ""
    seen: Dict[str, int] = {}
    in_table = False
    with open_soft(path, encoding) as handle:
        for number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if line.startswith("^"):
                if in_table:
                    study.warnings.append(f"Line {number}: table missing end marker")
                in_table = False
                match = ENTITY.fullmatch(line)
                if not match:
                    raise MetadataError(f"Line {number}: malformed entity header {line!r}")
                entity_kind, identifier = match.group(1).upper(), match.group(2)
                current, series_id = None, ""
                if entity_kind == "SERIES":
                    if not re.fullmatch(r"GSE\d+", identifier):
                        raise MetadataError(f"Line {number}: invalid GSE {identifier!r}")
                    series_id = identifier
                    if identifier not in study.series:
                        study.series.append(identifier)
                elif entity_kind == "SAMPLE":
                    if not re.fullmatch(r"GSM\d+", identifier):
                        raise MetadataError(f"Line {number}: invalid GSM {identifier!r}")
                    if identifier in seen:
                        raise MetadataError(
                            f"Duplicate GSM {identifier}: SAMPLE sections at lines "
                            f"{seen[identifier]} and {number}; refusing ambiguous merge")
                    seen[identifier] = number
                    current = Sample(identifier, number)
                    study.samples.append(current)
                continue
            marker = TABLE.fullmatch(line)
            if marker:
                beginning = marker.group(1).lower() == "begin"
                if beginning == in_table:
                    study.warnings.append(f"Line {number}: unexpected table marker")
                in_table = beginning
                continue
            if in_table or not line.strip():
                continue
            if line.startswith("!") and "=" in line:
                key, value = line[1:].split("=", 1)
                key = key.strip()
                value = value[1:] if value.startswith(" ") else value
                lower_key = key.lower()
                if current is not None and lower_key.startswith("sample_"):
                    current.metadata.append((key, value))
                elif entity_kind == "SERIES" and lower_key == "series_sample_id":
                    gsm = value.strip()
                    if re.fullmatch(r"GSM\d+", gsm):
                        study.membership[gsm].append(series_id)
                    else:
                        study.warnings.append(f"Line {number}: invalid Series_sample_id {value!r}")
                elif current is not None:
                    current.unparsed.append((number, line))
                    study.warnings.append(f"Line {number}: unexpected field in {current.gsm}")
            elif current is not None and not line.startswith("#"):
                current.unparsed.append((number, line))
                study.warnings.append(f"Line {number}: unparsed metadata in {current.gsm}")
    if in_table:
        study.warnings.append("End of file: table missing end marker")
    if not study.samples:
        raise MetadataError("No valid ^SAMPLE sections were found")
    return study


def split_characteristic(value: str) -> Optional[Tuple[str, str]]:
    """Accept Tag: Value, keeping further colons in Value."""
    label, separator, content = value.partition(":")
    if not separator or not label.strip():
        return None
    return label.strip(), content.strip()


def allocate_columns(names: Sequence[str], prefix: str) -> Dict[str, str]:
    """Disambiguate distinct original labels with deterministic suffixes."""
    groups: Dict[str, List[str]] = defaultdict(list)
    for name in sorted(set(names)):
        groups[slug(name)].append(name)
    result: Dict[str, str] = {}
    used = set()
    for base, originals in sorted(groups.items()):
        for original in originals:
            candidate = prefix + base
            if len(originals) > 1 or candidate in used:
                candidate += "__" + hashlib.sha256(original.encode("utf-8")).hexdigest()
            while candidate in used:
                candidate += "_"
            used.add(candidate)
            result[original] = candidate
    return result


def make_rows(study: Study) -> Tuple[List[Dict[str, str]], List[str]]:
    """Build a union schema without mixing series/platform fields into samples."""
    raw_names = [key for sample in study.samples for key, _ in sample.metadata]
    raw_columns = allocate_columns(raw_names, "geo__")
    labels: Dict[str, List[str]] = defaultdict(list)
    for sample in study.samples:
        for key, value in sample.metadata:
            match = CHARACTERISTIC.fullmatch(key.lower()[7:])
            parsed = split_characteristic(value) if match else None
            if match and parsed:
                labels[match.group(1)].append(parsed[0])
    char_columns = {channel: allocate_columns(names, f"characteristic_ch{channel}__")
                    for channel, names in labels.items()}
    rows: List[Dict[str, str]] = []
    for sample in study.samples:
        values: Dict[str, List[str]] = defaultdict(list)
        provenance: Dict[str, List[int]] = defaultdict(list)
        original: Dict[str, List[str]] = defaultdict(list)
        dynamic: Dict[str, List[str]] = defaultdict(list)
        explicit_series: List[str] = []

        def add(column: str, value: str, index: int) -> None:
            values[column].append(value.strip())
            provenance[column].append(index)

        for index, (key, raw_value) in enumerate(sample.metadata):
            original[key].append(raw_value)
            base = key.lower()[7:]  # Parser admits Sample_* fields only.
            value = raw_value.strip()
            for column, source in DIRECT_FIELDS.items():
                if base == source or (not re.search(r"_ch\d+$", source)
                                       and base == source + "_ch1"):
                    add(column, raw_value, index)
            if base == "source_name_ch1":
                add("tissue_source", raw_value, index)
            if re.fullmatch(r"supplementary_file(?:_\d+)?", base):
                add("supplementary_files", raw_value, index)
            if base == "geo_accession" and value != sample.gsm:
                raise MetadataError(f"{sample.gsm}: conflicting Sample_geo_accession {value!r}")
            if base == "series_id":
                if not re.fullmatch(r"GSE\d+", value):
                    study.warnings.append(f"{sample.gsm}: invalid Sample_series_id {value!r}")
                else:
                    explicit_series.append(value)
            if base in {"relation", "biosample", "sra", "biosample_accession",
                        "sra_experiment", "sra_experiment_accession",
                        "sra_run", "sra_project", "bioproject", "bioproject_accession"}:
                for column, pattern in ACCESSIONS.items():
                    for accession in pattern.findall(raw_value):
                        add(column, accession, index)
            match = CHARACTERISTIC.fullmatch(base)
            if match:
                parsed = split_characteristic(raw_value)
                if parsed:
                    label, content = parsed
                    dynamic[char_columns[match.group(1)][label]].append(content)
                    if match.group(1) == "1":
                        for column, aliases in ALIASES.items():
                            if slug(label) in aliases:
                                add(column, content, index)
                else:
                    dynamic[f"geo_unparsed_characteristics_ch{match.group(1)}_json"].append(raw_value)
                    study.warnings.append(f"{sample.gsm}: characteristic lacks a nonempty 'Tag:'")
        linked = study.membership.get(sample.gsm, [])
        associations = list(dict.fromkeys(linked + explicit_series))
        sources = []
        if linked:
            sources.append("Series_sample_id")
        if explicit_series:
            sources.append("Sample_series_id")
        if not associations and len(study.series) == 1:
            associations = study.series[:]
            sources.append("sole_family_series")
        if not associations:
            study.warnings.append(f"{sample.gsm}: no unambiguous GSE association")
        row = {name: cell(values.get(name, [])) for name in STANDARD_FIELDS}
        row.update(GSE=cell(associations), GSM=sample.gsm)
        row.update({raw_columns[key]: to_json(items) for key, items in original.items()})
        row.update({key: (to_json(items) if key.startswith("geo_unparsed_characteristics_") else cell(items))
                    for key, items in dynamic.items()})
        row["geo_metadata_json"] = to_json(sample.metadata)
        row["geo_unparsed_lines_json"] = to_json(sample.unparsed)
        row["standardized_provenance_json"] = to_json(provenance)
        row["gse_association_source"] = cell(sources)
        rows.append(row)
    extras = sorted({key for row in rows for key in row} - set(STANDARD_FIELDS))
    return rows, list(STANDARD_FIELDS) + extras


def validate(study: Study, rows: List[Dict[str, str]], expected: Optional[int]) -> Dict[str, int]:
    """Check identity, membership, expected count, and important-field coverage."""
    identifiers = [row["GSM"] for row in rows]
    if len(identifiers) != len(set(identifiers)):
        raise MetadataError("Internal validation failed: duplicate output GSM identifiers")
    if expected is None and len(study.series) == 1:
        expected = KNOWN_COUNTS.get(study.series[0])
    if expected is not None and len(rows) != expected:
        raise MetadataError(f"Expected {expected} samples; extracted {len(rows)}")
    missing_members = sorted(set(study.membership) - set(identifiers))
    if missing_members:
        study.warnings.append("Series references absent SAMPLE sections: " + ", ".join(missing_members))
    missing = {name: 0 for name in STANDARD_FIELDS}
    # Coverage uses underlying values, so repeated empty strings count as missing.
    for sample, row in zip(study.samples, rows):
        provenance = json.loads(row["standardized_provenance_json"])
        absent = set()
        for name in STANDARD_FIELDS:
            populated = bool(row[name].strip())
            if name in provenance:
                text = row[name]
                items = json.loads(text) if len(provenance[name]) > 1 else [text]
                populated = any(item.strip() for item in items)
            missing[name] += int(not populated)
            if not populated:
                absent.add(name)
        for name in IMPORTANT_FIELDS:
            if name in absent:
                study.warnings.append(f"{sample.gsm}: missing important field {name}")
    return missing


def write_csv(path: Path, rows: List[Dict[str, str]], columns: List[str]) -> None:
    """Atomically publish a UTF-8 CSV; do not leave a partially written output."""
    temporary: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=path.parent, prefix=f".{path.name}.",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(path))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def write_json(path: Path, payload: object) -> None:
    """Write QC report atomically without modifying another output on failure."""
    temporary: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=path.parent, prefix=f".{path.name}.",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(path))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def write_raw_jsonl(path: Path, study: Study, rows: List[Dict[str, str]]) -> None:
    """Preserve original sample fields in source order, one sample per line."""
    temporary: Optional[Path] = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=path.parent, prefix=f".{path.name}.",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            for sample, row in zip(study.samples, rows):
                record = {
                    "GSE": row["GSE"], "GSM": sample.gsm,
                    "sample_section_line": sample.line,
                    "sample_metadata": sample.metadata,
                    "unparsed_lines": sample.unparsed,
                    "gse_association_source": row["gse_association_source"],
                    "standardized_provenance": json.loads(row["standardized_provenance_json"]),
                }
                handle.write(to_json(record) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(path))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True, type=Path, help="Input .soft or .soft.gz")
    parser.add_argument("--output", required=True, type=Path, help="Output UTF-8 CSV")
    parser.add_argument("--expected-samples", type=int, help="Required total GSM count")
    parser.add_argument("--encoding", default="utf-8-sig", help="Input encoding (default: utf-8-sig)")
    parser.add_argument("--strict", action="store_true", help="Fail on any validation warning")
    parser.add_argument("--list-fields", action="store_true", help="Print every available CSV column")
    parser.add_argument("--qc-report", type=Path, help="Write validation report to JSON")
    parser.add_argument("--raw-jsonl", type=Path, help="Write original sample metadata to JSONL")
    args = parser.parse_args(argv)
    if args.expected_samples is not None and args.expected_samples < 1:
        parser.error("--expected-samples must be positive")
    try:
        source, target = args.input.expanduser().resolve(), args.output.expanduser().resolve()
        qc_path = args.qc_report.expanduser().resolve() if args.qc_report else None
        raw_path = args.raw_jsonl.expanduser().resolve() if args.raw_jsonl else None
        destinations = [target] + [p for p in (qc_path, raw_path) if p is not None]
        if source in destinations or len(set(destinations)) != len(destinations):
            raise MetadataError("Input, CSV, QC report, and JSONL paths must all be different")
        for output_path in destinations:
            if output_path.exists() and os.path.samefile(source, output_path):
                raise MetadataError("Output must not overwrite the input")
        study = parse_soft(source, args.encoding)
        rows, columns = make_rows(study)
        missing = validate(study, rows, args.expected_samples)
        for message in study.warnings[:20]:
            print(f"WARNING: {message}", file=sys.stderr)
        if len(study.warnings) > 20:
            print(f"WARNING: {len(study.warnings) - 20} additional findings", file=sys.stderr)
        if args.strict and study.warnings:
            raise MetadataError(f"Strict validation failed: {len(study.warnings)} warning(s)")
        for output_path in destinations:
            output_path.parent.mkdir(parents=True, exist_ok=True)
        # All parsing and strict validation occurs before any output is published.
        # Separate outputs are individually atomic, not a multi-file transaction.
        if raw_path is not None:
            write_raw_jsonl(raw_path, study, rows)
        if qc_path is not None:
            report = {
                "input": str(source), "GSE": study.series,
                "sample_count": len(rows), "expected_sample_count": (
                    args.expected_samples if args.expected_samples is not None else
                    KNOWN_COUNTS.get(study.series[0]) if len(study.series) == 1 else None),
                "duplicate_gsms": 0, "metadata_columns": columns,
                "missing_standardized_fields": missing,
                "warnings": study.warnings,
            }
            write_json(qc_path, report)
        write_csv(target, rows, columns)
        available = [name for name in STANDARD_FIELDS if missing[name] < len(rows)]
        print(f"Samples extracted: {len(rows)}; duplicate GSMs: 0")
        print("GSE records: " + (", ".join(study.series) or "none"))
        print(f"Metadata columns: {len(columns)}; validation warnings: {len(study.warnings)}")
        print("Available standardized fields: " + ", ".join(available))
        missing_text = ", ".join(f"{name}={count}/{len(rows)}"
                                 for name, count in missing.items() if count)
        print("Missing standardized values: " + (missing_text or "none"))
        if args.list_fields:
            print("All columns: " + ", ".join(columns))
        print(f"Output: {target}")
        if qc_path is not None:
            print(f"QC report: {qc_path}")
        if raw_path is not None:
            print(f"Raw metadata: {raw_path}")
        return 0
    except (MetadataError, OSError, EOFError, UnicodeError, LookupError, csv.Error,
            zlib.error) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
