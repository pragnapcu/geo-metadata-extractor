#!/usr/bin/env python3
"""Map GEO GSMs to author-supplied XLSX count columns (stdlib-only, Python 3.8+).

Reads an XLSX worksheet header, optional reviewed aliases, and GEO Sample_description retained
in the merged metadata CSV (typically geo__sample_description, a JSON array).
Matches identifiers exactly, preserving order and names. Does not process counts.
"""
import argparse
import csv
import json
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from posixpath import normpath

M = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P = 'http://schemas.openxmlformats.org/package/2006/relationships'
Q = lambda tag: f'{{{M}}}{tag}'


def col_number(ref):
    letters = re.match(r'^[A-Z]+', (ref or '').upper())
    if not letters:
        raise ValueError(f'Invalid Excel cell reference: {ref!r}')
    number = 0
    for letter in letters.group():
        number = number * 26 + ord(letter) - 64
    return number


def shared_strings(z):
    if 'xl/sharedStrings.xml' not in z.namelist():
        return []
    root = ET.fromstring(z.read('xl/sharedStrings.xml'))
    return [''.join(t.text or '' for t in item.iter(Q('t')))
            for item in root.findall(Q('si'))]


def worksheet_path(z, requested):
    book = ET.fromstring(z.read('xl/workbook.xml'))
    relationships = ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))
    targets = {item.get('Id'): item.get('Target')
               for item in relationships.findall(f'{{{P}}}Relationship')}
    sheets = [(item.get('name'), item.get(f'{{{R}}}id'))
              for item in book.find(Q('sheets'))]
    if not sheets:
        raise ValueError('Workbook contains no worksheets')
    if requested is None:
        requested = sheets[0][0]
    ids = dict(sheets)
    if requested not in ids:
        raise ValueError(f'Worksheet {requested!r} not found. Available: {list(ids)}')
    target = targets[ids[requested]]
    path = (target.lstrip('/') if target.startswith('/')
            else normpath('xl/' + target))
    if path not in z.namelist():
        raise ValueError(f'Worksheet XML missing: {path}')
    return requested, path


def cell_value(c, strings):
    kind = c.get('t')
    if kind == 'inlineStr':
        return ''.join(t.text or '' for t in c.iter(Q('t')))
    v = c.find(Q('v'))
    if v is None:
        return ''
    if kind == 's':
        idx = int(v.text)
        if idx < 0 or idx >= len(strings):
            raise ValueError(f'Invalid shared-string index {idx} at {c.get("r")}')
        return strings[idx]
    return v.text or ''


def get_header(path, sheet, header_row):
    with zipfile.ZipFile(path) as z:
        strings = shared_strings(z)
        name, member = worksheet_path(z, sheet)
        with z.open(member) as handle:
            for _, elem in ET.iterparse(handle, events=('end',)):
                if elem.tag != Q('row'):
                    continue
                rownum = int(elem.get('r', '0'))
                if rownum == header_row:
                    cells = {}
                    for c in elem.findall(Q('c')):
                        n = col_number(c.get('r'))
                        cells[n] = cell_value(c, strings)
                    return name, cells
                elem.clear()
    raise ValueError(f'Header row {header_row} not found')


def read_description(value):
    """CSV exported by extractor stores original SOFT values as JSON arrays."""
    if not value:
        return []
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError:
        return [value]
    if isinstance(parsed, list):
        return [item for item in parsed if isinstance(item, str) and item]
    if isinstance(parsed, str):
        return [parsed] if parsed else []
    raise ValueError(f'Unexpected description representation: {value!r}')


def read_aliases(path):
    """Read explicitly reviewed exceptions. Never silently normalize identifiers."""
    if path is None:
        return {}
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        required = {'GSM', 'count_column', 'review_status', 'evidence'}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f'Alias CSV requires columns: {sorted(required)}')
        aliases = {}
        for lineno, row in enumerate(reader, 2):
            gsm = row['GSM'].strip()
            name = row['count_column'].strip()
            status = row['review_status'].strip().lower()
            evidence = row['evidence'].strip()
            if not re.fullmatch(r'GSM\d+', gsm) or not name:
                raise ValueError(f'Alias row {lineno}: invalid GSM or empty count_column')
            if gsm in aliases:
                raise ValueError(f'Alias row {lineno}: duplicate GSM {gsm}')
            if status not in {'verified', 'assumed'} or not evidence:
                raise ValueError(f'Alias row {lineno}: {gsm} requires review_status=verified or assumed and nonempty evidence')
            aliases[gsm] = (name, status, evidence)
    return aliases


def build_mapping(metadata, header, description_field, aliases=None):
    aliases = aliases or {}
    if not header or min(header) != 1:
        raise ValueError('Expected gene identifier header in Excel column A')
    countcols = {number: val for number, val in header.items()
                 if number > 1 and val.strip()}
    if not countcols:
        raise ValueError('No sample headers found after column A')
    names = Counter(countcols.values())
    duplicates = sorted(name for name, total in names.items() if total > 1)
    by_name = defaultdict(list)
    for number, value in countcols.items():
        by_name[value].append(number)

    results, problems, matchedcols = [], [], set()
    known_gsms = {row.get('GSM', '').strip() for row in metadata}
    for gsm in sorted(set(aliases) - known_gsms):
        problems.append(f'Alias for unknown GSM: {gsm}')
    for row in metadata:
        gsm = row.get('GSM', '').strip()
        descriptions = read_description(row.get(description_field, ''))
        unique = list(dict.fromkeys(descriptions))
        hits = [(name, number) for name in unique for number in by_name.get(name, [])]
        chosen = None
        method, status, evidence = '', '', ''
        if len(hits) == 1:
            chosen = hits[0]
            method, status = 'exact_description', 'matched'
            if gsm in aliases and aliases[gsm][0] != chosen[0]:
                problems.append(f'{gsm}: alias conflicts with exact description match')
        elif len(hits) > 1:
            problems.append(f'{gsm}: ambiguous exact descriptions; matches={hits!r}')
        elif gsm in aliases:
            alias, status, evidence = aliases[gsm]
            alias_hits = by_name.get(alias, [])
            if len(alias_hits) == 1:
                chosen = (alias, alias_hits[0])
                method = 'reviewed_alias' if status == 'verified' else 'assumed_alias'
                status = 'matched' if method == 'reviewed_alias' else 'matched_assumed'
            else:
                problems.append(f'{gsm}: alias {alias!r} found {len(alias_hits)} times in count header')
        else:
            problems.append(f'{gsm}: unmatched; descriptions={unique!r}')
        if chosen is not None:
            name, number = chosen
            if number in matchedcols:
                problems.append(f'{gsm}: count column {number} assigned to more than one GSM')
            matchedcols.add(number)
        else:
            name, number = '', ''
            status = 'unmatched'
        results.append({
            'GSE': row.get('GSE', ''), 'GSM': gsm,
            'sample_title': row.get('sample_title', ''),
            'original_descriptions_json': json.dumps(descriptions, ensure_ascii=False),
            'internal_sample_id': name,
            'excel_column_number': number,
            'disease_state': row.get('disease_state', ''),
            'stimulation': row.get('stimulation', ''),
            'tissue_source': row.get('tissue_source', ''),
            'SRX': row.get('SRX', ''), 'SRR': row.get('SRR', ''),
            'mapping_method': method,
            'review_status': (aliases[gsm][1] if method in {'reviewed_alias', 'assumed_alias'} else 'not_required'),
            'evidence': evidence,
            'mapping_status': status,
        })
    for value in duplicates:
        problems.append(f'Duplicate count header: {value!r}')
    for number in sorted(set(countcols) - matchedcols):
        problems.append(f'Unmatched count column {number}: {countcols[number]!r}')
    return results, problems, len(countcols)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--metadata', required=True, type=Path, help='Merged GEO-SRA metadata CSV')
    ap.add_argument('--counts', required=True, type=Path, help='Author-supplied XLSX')
    ap.add_argument('--output', required=True, type=Path, help='Output mapping CSV')
    ap.add_argument('--sheet', default=None, help='Worksheet (default first)')
    ap.add_argument('--header-row', type=int, default=1)
    ap.add_argument('--description-field', default='geo__sample_description')
    ap.add_argument('--aliases', type=Path, help='Alias CSV with verified or assumed status (GSM,count_column,review_status,evidence)')
    args = ap.parse_args()
    try:
        if args.header_row < 1:
            raise ValueError('--header-row must be at least 1')
        with args.metadata.open(newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            if not reader.fieldnames or 'GSM' not in reader.fieldnames:
                raise ValueError('Metadata CSV needs a GSM column')
            if args.description_field not in reader.fieldnames:
                raise ValueError(f'{args.description_field!r} not found. Available description fields: '
                                 f'{[n for n in reader.fieldnames if "description" in n.lower()]}')
            metadata = list(reader)
        if not metadata:
            raise ValueError('Metadata CSV contains no rows')
        gsmcounts = Counter(row['GSM'] for row in metadata)
        repeated = [gsm for gsm, n in gsmcounts.items() if n > 1]
        if repeated:
            raise ValueError(f'Metadata has multiple rows per GSM: {repeated}. '
                             'Use sample-level metadata or explicitly collapse multi-run records.')
        sheet, header = get_header(args.counts, args.sheet, args.header_row)
        aliases = read_aliases(args.aliases)
        mapping, problems, num_countcols = build_mapping(metadata, header, args.description_field, aliases)
        print(f'Worksheet: {sheet}; GEO samples: {len(metadata)}; count sample columns: {num_countcols}')
        print(f'Exact description matches: {sum(r["mapping_method"] == "exact_description" for r in mapping)}')
        print(f'Reviewed alias matches: {sum(r["mapping_method"] == "reviewed_alias" for r in mapping)}')
        print(f'Assumed alias matches: {sum(r["mapping_method"] == "assumed_alias" for r in mapping)}')
        if problems:
            print('FAIL: Mapping was not one-to-one. No mapping CSV written.', file=sys.stderr)
            for message in problems:
                print(' - ' + message, file=sys.stderr)
            return 1
        columns = list(mapping[0])
        args.output.parent.mkdir(parents=True, exist_ok=True)
        if args.output.resolve() in (args.counts.resolve(), args.metadata.resolve()):
            raise ValueError('Output path must differ from input paths')
        with args.output.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(mapping)
        print(f'PASS: All {len(metadata)} samples mapped one-to-one; output: {args.output}')
        return 0
    except (OSError, ValueError, KeyError, zipfile.BadZipFile, ET.ParseError) as error:
        print(f'ERROR: {error}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
