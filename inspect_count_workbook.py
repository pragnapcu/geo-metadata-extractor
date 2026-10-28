
#!/usr/bin/env python3

import sys
import zipfile
import xml.etree.ElementTree as ET

path = sys.argv[1]

MAIN = (
    "http://schemas.openxmlformats.org/"
    "spreadsheetml/2006/main"
)
REL = (
    "http://schemas.openxmlformats.org/"
    "officeDocument/2006/relationships"
)
PKG = (
    "http://schemas.openxmlformats.org/"
    "package/2006/relationships"
)

with zipfile.ZipFile(path) as z:
    strings = []

    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for item in root.findall(f"{{{MAIN}}}si"):
            strings.append("".join(
                t.text or ""
                for t in item.iter(f"{{{MAIN}}}t")
            ))

    workbook = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(
        z.read("xl/_rels/workbook.xml.rels")
    )

    targets = {
        r.get("Id"): r.get("Target")
        for r in rels.findall(f"{{{PKG}}}Relationship")
    }

    sheets = workbook.find(f"{{{MAIN}}}sheets")

    for sheet in sheets:
        name = sheet.get("name")
        rid = sheet.get(f"{{{REL}}}id")
        target = targets[rid]

        if target.startswith("/"):
            member = target.lstrip("/")
        elif target.startswith("xl/"):
            member = target
        else:
            member = "xl/" + target

        print(f"\nWorksheet: {name}")

        with z.open(member) as stream:
            count = 0

            for event, elem in ET.iterparse(
                stream, events=("end",)
            ):
                if elem.tag != f"{{{MAIN}}}row":
                    continue

                values = []

                for c in elem.findall(f"{{{MAIN}}}c"):
                    v = c.find(f"{{{MAIN}}}v")

                    if c.get("t") == "inlineStr":
                        value = "".join(
                            t.text or ""
                            for t in c.iter(f"{{{MAIN}}}t")
                        )
                    elif v is None:
                        value = ""
                    elif c.get("t") == "s":
                        value = strings[int(v.text)]
                    else:
                        value = v.text or ""

                    values.append(
                        f"{c.get('r')}={value[:70]}"
                    )

                print(" | ".join(values[:18]))
                count += 1
                elem.clear()

                if count >= 5:
                    break
