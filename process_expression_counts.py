#!/usr/bin/env python3
"""Validate a GEO raw-count XLSX workbook and export a GSM-labelled TSV + QC.

Python 3.8+ standard library only. Reads an existing, one-to-one mapping CSV
from map_expression_columns_with_assumptions.py. No normalization is performed.
Assumed aliases remain explicitly flagged in the QC outputs.
"""
import argparse
import csv
import json
import math
import os
import re
import sys
import tempfile
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from posixpath import normpath

M = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P = 'http://schemas.openxmlformats.org/package/2006/relationships'
Q = lambda name: f'{{{M}}}{name}'


def excel_column(ref):
    m = re.match(r'^([A-Za-z]+)[1-9][0-9]*$', ref or '')
    if not m:
        raise ValueError(f'Invalid Excel reference: {ref!r}')
    n = 0
    for ch in m.group(1).upper():
        n = n * 26 + ord(ch) - 64
    return n


def read_mapping(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        required = {'GSM', 'internal_sample_id', 'excel_column_number', 'mapping_method'}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f'Mapping CSV must contain {sorted(required)}')
        rows = list(reader)
    if not rows:
        raise ValueError('Mapping CSV is empty')
    gsms, cols = set(), set()
    for row in rows:
        gsm = row['GSM'].strip()
        if not re.fullmatch(r'GSM\d+', gsm) or gsm in gsms:
            raise ValueError(f'Invalid or duplicate GSM: {gsm!r}')
        gsms.add(gsm)
        try:
            col = int(row['excel_column_number'])
        except (TypeError, ValueError):
            raise ValueError(f'Invalid Excel column for {gsm}') from None
        if col < 2 or col in cols:
            raise ValueError(f'Duplicate or invalid Excel column for {gsm}: {col}')
        cols.add(col)
        if not row['internal_sample_id']:
            raise ValueError(f'Missing internal sample ID for {gsm}')
        if row.get('mapping_status', 'matched') != 'matched':
            raise ValueError(f'Non-matched mapping for {gsm}')
        if row['mapping_method'] not in {'exact_description', 'reviewed_alias', 'assumed_alias'}:
            raise ValueError(f'Unknown mapping method for {gsm}: {row["mapping_method"]!r}')
        if row['mapping_method'] == 'assumed_alias' and not row.get('evidence', '').strip():
            raise ValueError(f'Assumed alias for {gsm} must contain evidence')
        row['_col'] = col
    return rows


def worksheet(z, requested):
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    relations = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    targets = {x.get('Id'): x.get('Target') for x in relations.findall(f'{{{P}}}Relationship')}
    sheets = {s.get('name'): s.get(f'{{{R}}}id') for s in wb.find(Q('sheets'))}
    if not sheets:
        raise ValueError('Workbook contains no worksheets')
    name = requested or next(iter(sheets))
    if name not in sheets:
        raise ValueError(f'Worksheet {name!r} not found: {list(sheets)}')
    target = targets[sheets[name]]
    member = target.lstrip('/') if target.startswith('/') else normpath('xl/' + target)
    if member not in z.namelist():
        raise ValueError(f'Missing worksheet XML: {member}')
    return name, member


def shared_strings(z):
    if 'xl/sharedStrings.xml' not in z.namelist():
        return []
    # Shared strings can be large; store only text, not an XML tree.
    result = []
    with z.open('xl/sharedStrings.xml') as f:
        for _, element in ET.iterparse(f, events=('end',)):
            if element.tag == Q('si'):
                result.append(''.join(t.text or '' for t in element.iter(Q('t'))))
                element.clear()
    return result


def get_value(cell, strings):
    if cell.find(Q('f')) is not None:
        raise ValueError(f'Formula cell not supported at {cell.get("r")}; raw counts must be values')
    typ = cell.get('t')
    if typ == 'inlineStr':
        return ''.join(t.text or '' for t in cell.iter(Q('t')))
    val = cell.find(Q('v'))
    if val is None or val.text is None:
        return ''
    if typ == 's':
        idx = int(val.text)
        if idx < 0 or idx >= len(strings):
            raise ValueError(f'Invalid shared-string index at {cell.get("r")}')
        return strings[idx]
    if typ in {'e', 'b'}:
        raise ValueError(f'Unexpected Excel error/Boolean cell at {cell.get("r")}')
    return val.text


def iter_rows(z, member, strings):
    with z.open(member) as stream:
        for _, elem in ET.iterparse(stream, events=('end',)):
            if elem.tag != Q('row'):
                continue
            number = int(elem.get('r', '0'))
            values = {}
            for cell in elem.findall(Q('c')):
                col = excel_column(cell.get('r'))
                if col in values:
                    raise ValueError(f'Duplicate Excel column at row {number}: {col}')
                values[col] = get_value(cell, strings)
            yield number, values
            elem.clear()


def count_value(raw, rownum, col):
    if raw == '':
        raise ValueError(f'Missing count in row {rownum}, Excel column {col}')
    try:
        n = float(raw)
    except (TypeError, ValueError):
        raise ValueError(f'Non-numeric count at row {rownum}, column {col}: {raw!r}') from None
    if not math.isfinite(n) or n < 0 or not n.is_integer():
        raise ValueError(f'Invalid raw count at row {rownum}, column {col}: {raw!r}')
    return int(n)


def atomic_path(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', suffix='.tmp', dir=str(path.parent))
    os.close(fd)
    return Path(temporary)


def run(args):
    mapping = read_mapping(args.mapping)
    columns = sorted(row['_col'] for row in mapping)
    lookup = {row['_col']: row for row in mapping}
    order = [lookup[col] for col in columns]  # Output in original workbook order.
    samples = {row['GSM']: {'GSM': row['GSM'], 'original_column': row['internal_sample_id'],
                             'mapping_method': row['mapping_method'],
                             'review_status': row.get('review_status', ''),
                             'evidence': row.get('evidence', ''),
                             'total_counts': 0, 'detected_genes': 0, 'zero_count_genes': 0}
               for row in order}
    genes, gene_count, all_zero = set(), 0, 0
    temp = atomic_path(args.output)
    try:
        with zipfile.ZipFile(args.counts) as z:
            name, member = worksheet(z, args.sheet)
            strings = shared_strings(z)
            with temp.open('w', newline='', encoding='utf-8') as out:
                writer = csv.writer(out, delimiter='\t', lineterminator='\n')
                writer.writerow(['EnsemblID'] + [x['GSM'] for x in order])
                saw_header = False
                for rownum, cells in iter_rows(z, member, strings):
                    if rownum < args.header_row:
                        continue
                    if rownum == args.header_row:
                        if not cells.get(1, '').strip():
                            raise ValueError('Empty gene ID header in column A')
                        present = {col for col, value in cells.items() if col > 1 and value.strip()}
                        if present != set(columns):
                            raise ValueError(f'Workbook sample columns {sorted(present)} differ from mapping {columns}')
                        for col in columns:
                            if cells.get(col) != lookup[col]['internal_sample_id']:
                                raise ValueError(f'Header mismatch at Excel column {col}: '
                                                 f'{cells.get(col)!r} != {lookup[col]["internal_sample_id"]!r}')
                        saw_header = True
                        continue
                    if not saw_header:
                        raise ValueError(f'Header row {args.header_row} not found')
                    if not any(v.strip() for v in cells.values()):
                        continue  # Ignore completely empty rows.
                    unexpected = [c for c, v in cells.items() if c not in set(columns + [1]) and v.strip()]
                    if unexpected:
                        raise ValueError(f'Unexpected nonempty columns on row {rownum}: {unexpected}')
                    gene = cells.get(1, '').strip()
                    if not gene:
                        raise ValueError(f'Missing gene identifier at Excel row {rownum}')
                    if gene in genes:
                        raise ValueError(f'Duplicate gene identifier at row {rownum}: {gene}')
                    genes.add(gene)
                    values = [count_value(cells.get(col, ''), rownum, col) for col in columns]
                    gene_count += 1
                    if not any(values):
                        all_zero += 1
                    for item, value in zip(order, values):
                        stat = samples[item['GSM']]
                        stat['total_counts'] += value
                        stat['detected_genes'] += int(value > 0)
                        stat['zero_count_genes'] += int(value == 0)
                    writer.writerow([gene] + values)
        if gene_count == 0:
            raise ValueError('No gene rows found')
        if args.expected_samples is not None and len(order) != args.expected_samples:
            raise ValueError(f'Expected {args.expected_samples} samples, mapped {len(order)}')
        if len({r['GSM'] for r in order}) != len(order):
            raise ValueError('Duplicate output GSM columns')
        qc = {'study': args.gse, 'sheet': name, 'source_workbook': str(args.counts),
              'mapping_file': str(args.mapping), 'gene_rows': gene_count,
              'sample_columns': len(order), 'all_zero_genes': all_zero,
              'assumed_mappings': [dict(GSM=r['GSM'], count_column=r['internal_sample_id'],
                                        evidence=r.get('evidence', ''))
                                   for r in order if r['mapping_method'] == 'assumed_alias'],
              'sample_order': [r['GSM'] for r in order], 'count_type': 'author_submitted_raw_integer_counts',
              'normalization': 'none', 'notes': 'Assumed aliases are not independently verified.'}
        for target in (args.output, args.sample_qc, args.qc_report):
            if target.resolve() in (args.counts.resolve(), args.mapping.resolve()):
                raise ValueError('Output path must not overwrite inputs')
        if len({p.resolve() for p in (args.output, args.sample_qc, args.qc_report)}) != 3:
            raise ValueError('Output paths must be distinct')
        tmp_sample, tmp_json = atomic_path(args.sample_qc), atomic_path(args.qc_report)
        try:
            with tmp_sample.open('w', newline='', encoding='utf-8') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(next(iter(samples.values())).keys()))
                writer.writeheader()
                writer.writerows(samples[r['GSM']] for r in order)
            tmp_json.write_text(json.dumps(qc, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
            os.replace(temp, args.output)
            os.replace(tmp_sample, args.sample_qc)
            os.replace(tmp_json, args.qc_report)
        finally:
            for p in (tmp_sample, tmp_json):
                if p.exists():
                    p.unlink()
    finally:
        if temp.exists():
            temp.unlink()
    print(f'PASS: {gene_count} genes x {len(order)} GSM samples')
    print(f'All-zero genes: {all_zero}; assumed mappings: {len(qc["assumed_mappings"])}')
    print('Counts:', args.output)
    print('Sample QC:', args.sample_qc)
    print('QC JSON:', args.qc_report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', required=True, type=Path)
    parser.add_argument('--mapping', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--sample-qc', required=True, type=Path)
    parser.add_argument('--qc-report', required=True, type=Path)
    parser.add_argument('--sheet', default='rawCounts')
    parser.add_argument('--header-row', type=int, default=1)
    parser.add_argument('--expected-samples', type=int, default=14)
    parser.add_argument('--gse', default='GSE178994')
    args = parser.parse_args()
    try:
        if args.header_row < 1 or (args.expected_samples is not None and args.expected_samples < 1):
            raise ValueError('Row and expected sample counts must be positive')
        run(args)
        return 0
    except (ValueError, OSError, KeyError, zipfile.BadZipFile, ET.ParseError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
