"""Colour arithmetic and palette search for ClonoView.

The numeric core is vendored from ../../tools/validate_palette.py, itself a port of the
data-viz skill's validate_palette.js: same OKLab conversion, same Machado-Oliveira-
Fernandes (2009) severity-1.0 colour-vision transforms, same thresholds, delta-E as
Euclidean distance in OKLab times 100. It is copied rather than imported so the viewer
is self-contained, and `tests/test_colour.py` checks the copy still reproduces that
file's documented reference numbers.

What is added here is the part that was never kept: the *search*. The 11 profile colours
below were produced by a search whose code no longer exists, so extending the scheme to
a dataset's cell types -- Mathew has 16 -- had nothing to extend. `search()` reproduces
the method described in analysis/V_embedding_viz/01_plot_embedding.py: maximise the
worst pair over ALL pairs under simulated protanopia and deuteranopia, subject to the
lightness band and chroma floor, then spread hue with whatever freedom is left.

THE HONEST LIMIT. A scatter where any class can neighbour any other is the all-pairs
regime, and it does not scale. Eleven classes clear the >=8 CVD / >=15 normal-vision
gates; twelve do not. Sixteen is far past that, so Mathew's cell types CANNOT be made
safe by colour alone at any assignment, and no search will change that. What the search
still buys is the best achievable worst pair, and the measured value is recorded next to
the palette so the number is visible rather than implied. Everything downstream then
treats direct labels as mandatory and offers click-to-isolate, because those are the
encodings that actually carry identity at this count.
"""

import json
import math
import os
import threading

BAND = {"light": (0.43, 0.77), "dark": (0.48, 0.67)}   # OKLCH L
CHROMA_FLOOR = 0.10
CVD_TARGET, CVD_FLOOR = 8.0, 6.0
NORMAL_FLOOR = 15.0
CONTRAST_MIN = 3.0

MACHADO = {
    "protan": ((0.152286, 1.052583, -0.204868),
               (0.114503, 0.786281, 0.099216),
               (-0.003882, -0.048116, 1.051998)),
    "deutan": ((0.367322, 0.860646, -0.227968),
               (0.280085, 0.672501, 0.047413),
               (-0.011820, 0.042940, 0.968881)),
    "tritan": ((1.255528, -0.076749, -0.178779),
               (-0.078411, 0.930809, 0.147602),
               (0.004733, 0.691367, 0.303900)),
}

# ---- the ClonoView surface ---------------------------------------------------
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
CONTEXT = "#dcdbd4"      # cells and clones with no profile; background points

# MANAscore. Deliberately outside the categorical palettes: it is an overlay meaning
# "this cell/clone was independently called tumour-reactive", not another class, and it
# has to read on top of whatever colouring is already there. Both marks carry a second
# channel besides hue -- a white halo on cells, a thick ring on clones -- because the
# palettes already contain reds and hue alone would not separate them reliably.
MANA_RED = "#e00025"
GRID = "#e1e0d9"

# ---- profiles ----------------------------------------------------------------
# Carried verbatim from analysis/V_embedding_viz/01_plot_embedding.py, where they were
# selected and validated: worst-pair CVD dE 8.0 and normal-vision dE 15.3 over all 55
# pairs. Slot k colours profile k+1. Reused rather than re-searched so that a figure in
# the paper and the same profile in this viewer are the same colour.
PROFILE_COLOURS = ["#2a78d6", "#30d80c", "#902400", "#54c0d8", "#d83c6c", "#841878",
                   "#540ccc", "#cc7800", "#786c18", "#cc78fc", "#fc9090"]

# Single-hue blue ramp for genuinely continuous magnitude (assignment confidence).
SEQ_BLUE = ["#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]

# Red -> purple -> blue ramp for the ordered "which interval" variable on arrows. It
# carries hue AND monotone lightness rather than lightness alone, which a single-hue
# ramp could not do on thin overlapping strokes.
SEQ_TIME = ["#e8442c", "#d4308c", "#a032c4", "#6040cc", "#2a78d6", "#104281"]


# ---- sRGB / OKLab ------------------------------------------------------------
def hex2srgb(h):
    h = h.strip().lstrip("#")
    return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]


def s2lin(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def lin2s(c):
    c = min(1.0, max(0.0, c))
    return 12.92 * c if c <= 0.0031308 else 1.055 * c ** (1 / 2.4) - 0.055


def lin(h):
    return [s2lin(c) for c in hex2srgb(h)]


def rel_lum(h):
    r, g, b = lin(h)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def contrast(a, b):
    hi, lo = sorted((rel_lum(a), rel_lum(b)), reverse=True)
    return (hi + 0.05) / (lo + 0.05)


def oklab_from_lin(rgb):
    r, g, b = rgb
    l = (0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b) ** (1 / 3)
    m = (0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b) ** (1 / 3)
    s = (0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b) ** (1 / 3)
    return (0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s,
            1.9779984951 * l - 2.4285922050 * m + 0.4505937099 * s,
            0.0259040371 * l + 0.7827717662 * m - 0.8086757660 * s)


def oklch(h):
    L, a, b = oklab_from_lin(lin(h))
    return L, math.hypot(a, b), math.degrees(math.atan2(b, a)) % 360


def lch2hex(L, C, H):
    """OKLCH to hex, or None when the colour falls outside the sRGB gamut.

    The inverse of `oklab_from_lin`, so the search can walk a perceptual grid and only
    keep the points sRGB can actually show.
    """
    a, b = C * math.cos(math.radians(H)), C * math.sin(math.radians(H))
    l_ = L + 0.3963377774 * a + 0.2158037573 * b
    m_ = L - 0.1055613458 * a - 0.0638541728 * b
    s_ = L - 0.0894841775 * a - 1.2914855480 * b
    l, m, s = l_ ** 3, m_ ** 3, s_ ** 3
    r = 4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s
    g = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s
    bb = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s
    if min(r, g, bb) < -0.001 or max(r, g, bb) > 1.001:
        return None
    return "#%02x%02x%02x" % tuple(round(255 * lin2s(c)) for c in (r, g, bb))


def simulate(h, kind):
    r, g, b = lin(h)
    M = MACHADO[kind]
    return [min(1.0, max(0.0, M[i][0] * r + M[i][1] * g + M[i][2] * b)) for i in range(3)]


def delta_e(h1, h2, kind=None):
    a = oklab_from_lin(simulate(h1, kind) if kind else lin(h1))
    b = oklab_from_lin(simulate(h2, kind) if kind else lin(h2))
    return 100 * math.dist(a, b)


# ---- measurement -------------------------------------------------------------
def measure(pal, surface=SURFACE, pairs="all"):
    """Worst-pair distances and per-slot properties. `pairs="all"` is the regime that
    applies to a scatter, where any class can land beside any other."""
    idx = ([(i, i + 1) for i in range(len(pal) - 1)] if pairs == "adjacent"
           else [(i, j) for i in range(len(pal)) for j in range(i + 1, len(pal))])
    worst_cvd, worst_cvd_pair = math.inf, None
    worst_norm, worst_norm_pair = math.inf, None
    for i, j in idx:
        cvd = min(delta_e(pal[i], pal[j], "protan"), delta_e(pal[i], pal[j], "deutan"))
        nm = delta_e(pal[i], pal[j])
        if cvd < worst_cvd:
            worst_cvd, worst_cvd_pair = cvd, (i, j)
        if nm < worst_norm:
            worst_norm, worst_norm_pair = nm, (i, j)
    slots = []
    for c in pal:
        L, C, H = oklch(c)
        slots.append(dict(hex=c, L=L, C=C, H=H, contrast=contrast(c, surface)))
    return dict(n=len(pal), n_pairs=len(idx),
                worst_cvd=worst_cvd, worst_cvd_pair=worst_cvd_pair,
                worst_normal=worst_norm, worst_normal_pair=worst_norm_pair,
                low_contrast=[i for i, s in enumerate(slots)
                              if s["contrast"] < CONTRAST_MIN],
                slots=slots)


# ---- search ------------------------------------------------------------------
def _grid(mode="light", dL=0.055, dC=0.03, dH=7.5):
    lo, hi = BAND[mode]
    out = []
    L = lo
    while L <= hi + 1e-9:
        C = CHROMA_FLOOR
        while C <= 0.33:
            H = 0.0
            while H < 360:
                h = lch2hex(L, C, H)
                if h:
                    out.append(h)
                H += dH
            C += dC
        L += dL
    return sorted(set(out))


def search(n, mode="light", seed_hue=250.0, rounds=8, cvd_step=0.25,
           normal_floor=NORMAL_FLOOR):
    """Pick `n` colours by the method the profile palette was built with.

    The objective is lexicographic, and the ORDER of its terms is the whole design.
    Both were established by measurement in the original script and both are easy to
    get wrong:

      1. meet the normal-vision floor. Not a tiebreak -- a constraint. Maximising
         colour-vision distance alone happily returns an 11-slot palette whose worst
         normal-vision pair is dE 11.3, well under the floor of 15: it buys separation
         for colour-blind viewers by spending it for everyone else.
      2. maximise the worst pair under simulated protanopia and deuteranopia, over ALL
         pairs, quantised to `cvd_step`. The quantisation is what leaves any freedom
         for term 3; without it a continuous objective is never tied and hue never gets
         to act.
      3. spread hue with what remains. This has to come LAST. Ordering hue spread above
         the colour-vision term collapses the worst CVD pair from 8.8 to 2.8, because an
         even spread around the hue wheel lands slots directly on the protan/deutan
         confusion line. But omitting it entirely is just as bad in the other direction:
         with no hue term the 16-slot search returns nine near-blues that clear every
         numeric gate and are still hopeless to tell apart on a dense scatter.
      4. worst normal-vision pair, as a final tiebreak.

    Greedy farthest-point insertion, then repeated local replacement: each slot in turn
    is re-chosen against the other n-1 held fixed, keeping whatever improves the score.
    """
    cand = _grid(mode)
    sim = {c: (oklab_from_lin(simulate(c, "protan")),
               oklab_from_lin(simulate(c, "deutan")),
               oklab_from_lin(lin(c)), oklch(c)[2]) for c in cand}

    def d_cvd(a, b):
        pa, da, _, _ = sim[a]
        pb, db, _, _ = sim[b]
        return 100 * min(math.dist(pa, pb), math.dist(da, db))

    def d_norm(a, b):
        return 100 * math.dist(sim[a][2], sim[b][2])

    def hue_gap(pal):
        """Smallest circular gap between any two hues -- larger is better spread."""
        hs = sorted(sim[c][3] for c in pal)
        if len(hs) < 2:
            return 360.0
        gaps = [hs[i + 1] - hs[i] for i in range(len(hs) - 1)]
        gaps.append(360 - hs[-1] + hs[0])
        return min(gaps)

    def score(pal):
        pairs = [(i, j) for i in range(len(pal)) for j in range(i + 1, len(pal))]
        cv = min(d_cvd(pal[i], pal[j]) for i, j in pairs)
        nm = min(d_norm(pal[i], pal[j]) for i, j in pairs)
        if nm >= normal_floor:
            return (1, int(cv / cvd_step), hue_gap(pal), cv)
        # Below the floor the score climbs toward it instead. A plain binary feasible
        # flag stalls the local search: one slot moves at a time, and no single move
        # gets an infeasible palette over the floor while also holding CVD up, so the
        # search sits in the infeasible region reporting a good CVD number that is not
        # actually allowed. Measured: with a binary flag, n=11 converged to normal-
        # vision dE 11.3, worse than the shipped profile palette's 15.3.
        return (0, nm, 0.0, cv)

    def refine(pal):
        best = score(pal)
        for _ in range(rounds):
            moved = False
            for k in range(n):
                others = pal[:k] + pal[k + 1:]
                cur, cand_best, cand_score = pal[k], pal[k], best
                for c in cand:
                    if c in others:
                        continue
                    s = score(others + [c])
                    if s > cand_score:
                        cand_best, cand_score = c, s
                if cand_best != cur:
                    pal[k], best, moved = cand_best, cand_score, True
            if not moved:
                break
        return pal, best

    if n == 1:
        return [min(cand, key=lambda c: abs(oklch(c)[2] - seed_hue))]

    # Several starts, because the local search is not convex and a single greedy seed
    # lands in a different basin depending on where it begins.
    best_pal, best_score = None, None
    for h0 in (seed_hue, 30.0, 90.0, 150.0, 210.0, 330.0):
        pal = [min(cand, key=lambda c: abs(oklch(c)[2] - h0))]
        while len(pal) < n:
            pal.append(max(cand, key=lambda c: (min(d_cvd(c, p) for p in pal),
                                                min(d_norm(c, p) for p in pal))))
        pal, s = refine(pal)
        if best_score is None or s > best_score:
            best_pal, best_score = pal, s
    # Order by hue so neighbouring legend entries are neighbouring hues. The objective
    # is over all pairs, so ordering cannot change any measured value.
    return sorted(best_pal, key=lambda c: oklch(c)[2])


# ---- cell-type palettes ------------------------------------------------------
# Cell-type counts vary by dataset (Liu 3, Wang 9, Mathew 16), and a palette searched
# for one count is not a prefix of the palette for another -- the objective is over all
# pairs, so adding a class re-solves the whole problem. Results are therefore stored per
# count in palettes.json, alongside the measurement that was achieved, and regenerated
# by tools/build_palettes.py rather than being chosen by hand.
_PAL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "palettes.json")
_PAL_LOCK = threading.Lock()
_PAL = None


def _palettes():
    global _PAL
    if _PAL is None:
        if os.path.exists(_PAL_PATH):
            with open(_PAL_PATH) as fh:
                _PAL = json.load(fh)
        else:
            _PAL = {}
    return _PAL


def cell_type_palette(n):
    """`n` colours for a dataset's cell types.

    Served from palettes.json when it is there. If it is not -- a dataset with a cell
    -type count nobody has prepared before -- the search runs once and the result is
    written back, so the cost is paid on first use and never again. The search takes
    seconds for small n and about a minute for n around 16.
    """
    if n <= 0:
        return []
    with _PAL_LOCK:
        p = _palettes()
        key = str(n)
        if key not in p:
            pal = search(n)
            m = measure(pal)
            p[key] = dict(colours=pal, worst_cvd=round(m["worst_cvd"], 2),
                          worst_normal=round(m["worst_normal"], 2),
                          meets_cvd_target=m["worst_cvd"] >= CVD_TARGET,
                          meets_normal_floor=m["worst_normal"] >= NORMAL_FLOOR,
                          low_contrast=m["low_contrast"])
            tmp = _PAL_PATH + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(p, fh, indent=1, sort_keys=True)
            os.replace(tmp, _PAL_PATH)
        return list(p[key]["colours"])


def cell_type_palette_quality(n):
    """What the palette for `n` classes actually measures, so the UI can say so rather
    than implying a separation it does not have."""
    cell_type_palette(n)
    q = dict(_palettes()[str(n)])
    q.pop("colours", None)
    return q
