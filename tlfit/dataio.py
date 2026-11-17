"""Tolerant reading of transients (txt/dat/csv, decimal comma, headers)."""
from __future__ import annotations

import re

import numpy as np

from .i18n import Msg

_num = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def parse_table(text: str):
    """Returns (numeric matrix, column names). Ignores non-numeric lines."""
    lines = [ln.strip() for ln in text.replace("\r", "\n").split("\n")]
    lines = [ln for ln in lines if ln and not ln.startswith(("#", "%", "//"))]
    if not lines:
        raise ValueError(Msg("empty file"))
    sample = "\n".join(lines[:50])
    # separator: ';' or tab take priority; with a comma, decide whether it is decimal
    if ";" in sample:
        sep, decimal_comma = ";", "," in sample
    elif "," in sample:
        has_dot = "." in sample
        spaced = re.search(r"\d\s+[-+.\d]", sample) is not None
        if not has_dot and spaced:
            sep, decimal_comma = None, True        # "0,002   1,0121"
        else:
            sep, decimal_comma = ",", False        # "0.002,1.0121"
    else:
        sep, decimal_comma = None, False
    rows, header = [], None
    for ln in lines:
        parts = [p.strip() for p in (ln.split(sep) if sep else ln.split())]
        parts = [p for p in parts if p != ""]
        if decimal_comma:
            parts = [p.replace(",", ".") for p in parts]
        if parts and all(_num.match(p) for p in parts):
            rows.append([float(p) for p in parts])
        elif not rows and header is None:
            header = parts
    if not rows:
        raise ValueError(Msg("no numeric line found"))
    ncol = max(set(len(r) for r in rows), key=[len(r) for r in rows].count)
    data = np.array([r for r in rows if len(r) == ncol], float)
    if header is None or len(header) != ncol:
        header = [Msg("column {i}", i=i + 1) for i in range(ncol)]
    return data, header


# ---------------------------------------------------------------------------
# Excitation power of each transient (for series)
# ---------------------------------------------------------------------------

_UNIT = {"w": 1e3, "mw": 1.0, "uw": 1e-3, "µw": 1e-3}
_NUM = r"(\d+(?:[.,]\d+)?(?:[eE][+-]?\d+)?)"
_WITH_UNIT = re.compile(_NUM + r"\s*(mW|µW|uW|W)(?![a-zA-Z])", re.IGNORECASE)
_HEADER = re.compile(r"(?:pot[eê]ncia|power|\bpe?)\s*[:=]\s*" + _NUM + r"\s*(mW|µW|uW|W)?",
                     re.IGNORECASE)
_P_PREFIX = re.compile(r"(?:^|[^a-zA-Z])(?:p|pot|pe)[_\-=]?" + _NUM + r"(?![\d.])", re.IGNORECASE)


def _tofloat(s):
    return float(s.replace(",", "."))


def guess_power(name: str, text: str = ""):
    """Power in mW from the file name or from the header.

    Order: number with unit in the name ('glass_150mW.dat', '0,2W.txt'); header
    line ('# Power: 150 mW', '# Potência: 150 mW', 'P = 0.15 W'); number after 'P' in the
    name ('P150'). Without a unit, mW is assumed. Returns (value in mW or None, source),
    with source "name", "header" or "name_nounit" (translated by the interface)."""
    stem = re.sub(r"\.[A-Za-z0-9]+$", "", name)
    m = _WITH_UNIT.search(stem)
    if m:
        return _tofloat(m.group(1)) * _UNIT[m.group(2).lower()], "name"
    for ln in text.splitlines()[:30]:
        s = ln.strip()
        if not s or not (s[0] in "#%/" or re.match(r"[A-Za-zÀ-ú]", s)):
            continue
        h = _HEADER.search(s)
        if h:
            return _tofloat(h.group(1)) * _UNIT[(h.group(2) or "mW").lower()], "header"
    m = _P_PREFIX.search(stem)
    if m:
        return _tofloat(m.group(1)), "name_nounit"
    return None, ""
