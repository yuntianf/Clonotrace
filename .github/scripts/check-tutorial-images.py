"""Check that the published tutorial includes every bundled figure, unchanged."""

import hashlib
from html.parser import HTMLParser
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit


class Images(HTMLParser):
    def __init__(self):
        super().__init__()
        self.sources = []

    def handle_starttag(self, tag, attrs):
        if tag == "img":
            attrs = dict(attrs)
            # pkgdown may include an empty, decorative package-logo element.
            if not attrs.get("src") and "logo" in attrs.get("class", "").split():
                return
            self.sources.append(attrs.get("src", ""))


def check(repo, html):
    site = (repo / "docs").resolve()
    expected = set(re.findall(
        r'cached_demo_figure\("([^"\n]+)"\)',
        (repo / "vignettes/Clonotrace.Rmd").read_text(encoding="utf-8"),
    ))
    if not expected:
        raise ValueError("No bundled tutorial figures found in the source")
    images = Images()
    images.feed(html.read_text(encoding="utf-8"))
    found = set()
    for src in images.sources:
        url = urlsplit(src)
        if not src or url.scheme or url.netloc or url.path.startswith("/"):
            raise ValueError(f"Tutorial image is not a portable relative path: {src!r}")
        image = (html.parent / unquote(url.path)).resolve()
        if not image.is_relative_to(site) or not image.is_file():
            raise ValueError(f"Tutorial image is missing or outside the website: {src}")
        if image.name in expected:
            original = repo / "inst/extdata/demo_files" / image.name
            if hashlib.sha256(image.read_bytes()).digest() != hashlib.sha256(
                original.read_bytes()).digest():
                raise ValueError(f"Published figure differs from the original: {src}")
            found.add(image.name)
    if found != expected:
        raise ValueError(f"Missing tutorial figures: {sorted(expected - found)}")
    print(f"PASS: all {len(expected)} tutorial figures have portable links and "
          "are byte-identical to the bundled originals.")


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    page = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else (
        root / "docs/articles/Clonotrace.html")
    try:
        check(root, page)
    except (OSError, ValueError) as error:
        sys.exit(f"FAIL: {error}")
