#!/usr/bin/env python3
"""Check that every relative markdown link in the docs lands somewhere.

    .github/scripts/check_links.py

Only local links are followed. A link to another page is checked for the file,
and a link with a #fragment is checked against that page's headings too, since
a heading being renamed is the common way these break.

External URLs are left alone: reaching them needs the network, and a rate limit
or a slow host would fail this for reasons that have nothing to do with the
change under test.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAGES = sorted(ROOT.glob("*.md")) + sorted(ROOT.glob("docs/*.md")) + sorted(
    ROOT.glob("assets/**/*.md")
)

LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)\)")


def anchors(path):
    """GitHub's heading anchors: lowercase, punctuation dropped, spaces to dashes."""
    found = set()
    fenced = False
    for line in path.read_text().splitlines():
        if line.lstrip().startswith("```"):
            fenced = not fenced
        elif not fenced and line.startswith("#"):
            text = line.lstrip("#").strip().lower()
            found.add(re.sub(r"[^\w\s-]", "", text).replace(" ", "-"))
    return found


def main():
    problems = []
    for page in PAGES:
        for target in LINK.findall(page.read_text()):
            if target.startswith(("http://", "https://", "mailto:")):
                continue
            path, _, fragment = target.partition("#")
            if not path:
                if fragment not in anchors(page):
                    problems.append(f"{page.relative_to(ROOT)} -> #{fragment}")
                continue
            dest = (page.parent / path).resolve()
            if not dest.exists():
                problems.append(f"{page.relative_to(ROOT)} -> {target} (no such file)")
            elif fragment and dest.suffix == ".md" and fragment not in anchors(dest):
                problems.append(f"{page.relative_to(ROOT)} -> {target} (no such heading)")

    if problems:
        print("Broken links:")
        for p in problems:
            print(f"  {p}")
        sys.exit(1)
    print(f"{len(PAGES)} pages checked, all internal links resolve")


if __name__ == "__main__":
    main()
