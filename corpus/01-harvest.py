#!/usr/bin/env python3
"""
Tier 1 of the outcome corpus: version chains of Cochrane reviews, with the
"Authors' conclusions" of each version, harvested from Europe PMC's open REST
API.

Nothing here touches cochranelibrary.com. Everything below is served by
Europe PMC under its normal terms, which is what makes this the backbone of
the corpus rather than the part that needs anyone's permission.

The version chain does not need a lookup: a Cochrane DOI carries it. The
review is CD005343 and the version is the `.pubN` suffix, so grouping by
review id and sorting by N reconstructs the whole history from metadata that
is already in hand.

Writes corpus/data/versions.jsonl, one line per review version.
"""

import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
QUERY = 'JOURNAL:"Cochrane Database Syst Rev"'
PAGE = 1000
PAUSE = 0.4          # polite; Europe PMC asks for reasonable use
OUT = os.path.join(os.path.dirname(__file__), "data", "versions.jsonl")

# 10.1002/14651858.CD005343.pub7 -> ("CD005343", 7); an unsuffixed DOI is pub1.
DOI_RE = re.compile(r"10\.1002/14651858\.(CD\d+)(?:\.pub(\d+))?", re.I)
CONCL_RE = re.compile(r"<h4>Authors' conclusions</h4>(.*?)(?=<h4>|$)", re.S)


def fetch(cursor):
    params = {
        "query": QUERY,
        "resultType": "core",
        "format": "json",
        "pageSize": str(PAGE),
        "cursorMark": cursor,
    }
    url = BASE + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={
        "User-Agent": "staleness-corpus/0.1 (research; jnverbel@gmail.com)"
    })
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def strip_tags(s):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", s)).strip()


def parse(rec):
    doi = rec.get("doi") or ""
    m = DOI_RE.search(doi)
    if not m:
        return None
    review_id = m.group(1).upper()
    version = int(m.group(2)) if m.group(2) else 1

    abstract = rec.get("abstractText") or ""
    c = CONCL_RE.search(abstract)
    conclusions = strip_tags(c.group(1)) if c else None

    # A protocol is not a review: it has no results to have changed.
    ptypes = (rec.get("pubTypeList") or {}).get("pubType") or []
    is_protocol = any("protocol" in str(t).lower() for t in ptypes) or \
        "This is a protocol" in abstract[:400]

    return {
        "review_id": review_id,
        "version": version,
        "doi": doi,
        "pmid": rec.get("pmid"),
        "pmcid": rec.get("pmcid"),
        "year": rec.get("pubYear"),
        "date": rec.get("firstPublicationDate"),
        "title": rec.get("title"),
        "is_protocol": is_protocol,
        "has_conclusions": conclusions is not None,
        "conclusions": conclusions,
        "abstract": abstract,
    }


def main():
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    cursor, seen, pages = "*", 0, 0
    with open(OUT, "w") as fh:
        while True:
            d = fetch(cursor)
            hits = d.get("hitCount", 0)
            results = d.get("resultList", {}).get("result", [])
            if not results:
                break
            for rec in results:
                row = parse(rec)
                if row:
                    fh.write(json.dumps(row, ensure_ascii=False) + "\n")
                    seen += 1
            pages += 1
            nxt = d.get("nextCursorMark")
            print(f"  página {pages:3d}  ·  {seen:6d} versiones  de ~{hits}",
                  file=sys.stderr, flush=True)
            if not nxt or nxt == cursor:
                break
            cursor = nxt
            time.sleep(PAUSE)
    print(f"escritas {seen} versiones en {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
