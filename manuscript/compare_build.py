"""Compare every page, link target and figure in a new build with the released PDFs."""
import argparse
import hashlib
import json
from pathlib import Path

import pymupdf as pdf

HERE = Path(__file__).resolve().parent


def links(page):
    return [(link["kind"], link.get("uri"), link.get("page"), tuple(link["from"]),
             str(link.get("to")), link.get("nameddest")) for link in page.get_links()]


def compare(name, directory, dpi):
    released, rebuilt = HERE / name, directory / name
    with pdf.open(released) as left, pdf.open(rebuilt) as right:
        changed = {"text": [], "image": [], "links": []}
        nlinks = 0
        for index, (a, b) in enumerate(zip(left, right), 1):
            if a.get_text() != b.get_text():
                changed["text"].append(index)
            if a.rect != b.rect or a.get_pixmap(dpi=dpi).samples != b.get_pixmap(dpi=dpi).samples:
                changed["image"].append(index)
            alinks, blinks = links(a), links(b)
            nlinks += len(alinks)
            if alinks != blinks:
                changed["links"].append(index)
        return {"file": str(name), "released_pages": len(left), "rebuilt_pages": len(right),
                "link_annotations": nlinks, "dpi": dpi, "changed_pages": changed,
                "identical_pdf_bytes": hashlib.sha256(released.read_bytes()).digest() == hashlib.sha256(rebuilt.read_bytes()).digest(),
                "matches": len(left) == len(right) and not any(changed.values())}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=144)
    args = parser.parse_args()
    if args.dpi < 1:
        parser.error("--dpi must be positive")
    files = [Path("manuscript.pdf"), Path("supplement.pdf"), *sorted(p.relative_to(HERE) for p in (HERE / "figures").glob("*.pdf"))]
    report = [compare(name, args.build, args.dpi) for name in files]
    print(json.dumps(report, indent=2))
    if not all(item["matches"] for item in report):
        raise SystemExit("PDF pages differ. See file, comparison and page number above.")
