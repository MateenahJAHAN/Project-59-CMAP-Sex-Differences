"""Check repository-relative Markdown links without requesting third-party webpages."""
import argparse
from functools import cache
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r'\]\((<[^>]+>|[^\s)]+)(?:\s+["\'][^)"\']*["\'])?\)')
FENCES = re.compile(r'```.*?```', flags=re.DOTALL)


@cache
def markdown_anchors(path):
    markup = FENCES.sub("", path.read_text())
    ids = set(re.findall(r'<a\s+(?:id|name)=["\']([^"\']+)["\']', markup))
    used = set()
    for title in re.findall(r"^#{1,6}\s+(.+?)\s*$", markup, flags=re.MULTILINE):
        slug = re.sub(r"[^\w -]", "", title.lower()).replace(" ", "-")
        tag, suffix = slug, 1
        while tag in used:
            tag = f"{slug}-{suffix}"
            suffix += 1
        used.add(tag)
    return ids | used


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    errors, count, fragments = [], 0, 0
    for path in sorted(ROOT.rglob("*.md")):
        if ".git" in path.relative_to(ROOT).parts:
            continue
        markup = FENCES.sub("", path.read_text())
        for match in LINK.finditer(markup):
            link = urlsplit(match[1].strip("<>"))
            if link.scheme or link.netloc or not (link.path or link.fragment):
                continue
            count += 1
            target = path.parent / unquote(link.path) if link.path else path
            if not target.exists():
                errors.append(f"{path.relative_to(ROOT)} -> {link.path}")
            elif link.fragment and target.suffix.lower() == ".md":
                fragments += 1
                if unquote(link.fragment) not in markdown_anchors(target.resolve()):
                    errors.append(f"{path.relative_to(ROOT)} -> {link.path}#{link.fragment}")
    if errors:
        raise SystemExit("Missing local Markdown targets:\n" + "\n".join(errors))
    print(f"All {count} repository-relative Markdown links resolve, including {fragments} section anchors.")


if __name__ == "__main__":
    main()
