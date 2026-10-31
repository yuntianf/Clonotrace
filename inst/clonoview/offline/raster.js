/* The cell panel's backdrop, drawn in the browser instead of by matplotlib.
 *
 * This replaces render.py's `background()`, and it is the change that makes the offline
 * export possible at all. In the served version the backdrop is a PNG rendered on demand
 * and cached by a hash of (mode, isolate, profile filter, MANAscore, width). Those
 * combinations multiply: three colour modes, isolate any one of up to sixteen cell types
 * or eleven profiles, filter to any profile, MANAscore on or off, at any width the
 * browser window happens to be -- about 1,200 images per width for Mathew, at 0.7 MB
 * each. Pre-rendering them for an export is not a large job, it is an impossible one.
 * Drawing them here removes the question: nothing is enumerated, and the width is
 * whatever the window is.
 *
 * WHY ImageData AND NOT ctx.arc(). Wang has 273,508 cells. A path-and-fill per cell is a
 * canvas API call per cell and lands in the region of a second, which would be felt on
 * every change of filter -- the served version's own budget for this was 0.26-0.37 s of
 * matplotlib. Compositing a small pre-computed coverage stamp straight into an ImageData
 * buffer is about twenty pixel writes per cell with no API crossing, and finishes in tens
 * of milliseconds. The stamp is supersampled 4x4, so the marks are antialiased the way
 * matplotlib's are rather than being hard squares.
 *
 * The sizes, the ordering, the one-mask rule and the MANAscore halo are all render.py's,
 * and the reasoning for each is in that file. What is repeated here is only the parts a
 * reader of this file has to be able to check against it.
 */

window.CV = window.CV || {};

(function (CV) {
  'use strict';

  const DPI = 100.0;          // render.py's figure dpi; sets the point-to-pixel ratio
  const PT_PX = DPI / 72.0;   // 1 typographic point in pixels at that dpi

  /* render.py `_point_size`: area in points squared, scaled by (width/1400)^2 so a mark
   * covers the same fraction of the panel at any width. The floor exists because marks
   * below about 2 pt^2 are sub-pixel, and antialiasing then spreads each into a pale tint
   * of itself -- a fully saturated palette renders as pastel and the structure washes
   * out. */
  function pointSize(n, width) {
    const scale = Math.pow(width / 1400.0, 2);
    if (n > 200000) return 3.0 * scale;
    if (n > 100000) return 4.0 * scale;
    return 6.0 * scale;
  }

  /* matplotlib's scatter `s` is nominally the marker's area in points squared, so the
   * marker is sqrt(s) points across -- but only for marker 'o'. THE '.' MARKER IS HALF
   * THAT SIZE: markers.py gives it a `_point_size_reduction` of 0.5, so the same `s`
   * draws a dot of half the diameter and a QUARTER of the area. render.py uses '.' for
   * the cells and 'o' for the MANAscore overlay, so both factors are in play in one
   * picture.
   *
   * Missing this made every cell four times the area it should be. It did not look like
   * an error -- the map was the right shape, the right colours, in the right places, just
   * denser and more saturated -- and that is exactly why it is written down here rather
   * than left as a constant. */
  const POINT_REDUCTION = 0.5;                 // matplotlib markers.py, for '.'

  function radiusPx(s, marker) {
    const red = (marker === 'o') ? 1.0 : POINT_REDUCTION;
    return Math.sqrt(s) * 0.5 * red * PT_PX;
  }

  function hexToRgb(h) {
    const s = h.charAt(0) === '#' ? h.slice(1) : h;
    const v = s.length === 3
      ? [s[0] + s[0], s[1] + s[1], s[2] + s[2]]
      : [s.slice(0, 2), s.slice(2, 4), s.slice(4, 6)];
    return [parseInt(v[0], 16), parseInt(v[1], 16), parseInt(v[2], 16)];
  }

  /* A disc's per-pixel coverage, supersampled 4x4, for a disc sitting at a given
   * sub-pixel offset from the centre of the pixel it lands in.
   *
   * TWO THINGS HERE ARE EASY TO GET WRONG AND WERE, and both look like a rendering
   * difference rather than a bug:
   *
   *   The supersample offsets must be centred on the pixel, spanning -0.5 to +0.5. An
   *   earlier version spanned 0 to +1, which shifts every mark half a pixel down and
   *   right AND halves the coverage of the pixel a mark is centred in -- so each mark
   *   spread over four pixels at a quarter alpha instead of landing on one. The picture
   *   was the right shape with the right colours and simply looked softer; the measured
   *   symptom was 46% ink coverage against matplotlib's 31% while the saturated pixels of
   *   most classes fell to nearly zero.
   *
   *   The sub-pixel offset has to be carried, not rounded away. Snapping every cell to
   *   the nearest pixel centre quantises 273,508 positions onto a grid and visibly
   *   coarsens dense regions. Agg places marks continuously; quantising the offset to
   *   quarters of a pixel gets within a rounding error of that for sixteen cached stamps
   *   per radius.
   */
  /* 8 sub-pixel positions per axis, not 4. The bucket is what stands in for a mark's true
   * position, so its width is a systematic error: at 4 buckets a mark is never modelled
   * as sitting across a pixel boundary, which is exactly the position that produces two
   * half-lit pixels and no saturated one. Sparse classes -- the ones drawn as a few
   * thousand non-overlapping marks -- came out with 30-60% more saturated pixels than
   * matplotlib for that reason alone. Stamps are cached per (radius, bucket), so 64 of
   * them per radius costs a few hundred microseconds once. */
  const SUBPIX = 8;
  const stampCache = new Map();
  function stamp(r, bx, by) {
    const key = Math.round(r * 64) * (SUBPIX * SUBPIX) + bx * SUBPIX + by;
    const hit0 = stampCache.get(key);
    if (hit0) return hit0;
    const rr = Math.max(r, 0.35);
    // The offset this bucket represents, in pixels from the pixel's centre.
    const fx = -0.5 + (bx + 0.5) / SUBPIX;
    const fy = -0.5 + (by + 0.5) / SUBPIX;
    const half = Math.ceil(rr + 1.0);
    const size = 2 * half + 1;
    const a = new Float32Array(size * size);
    /* 16x16 samples per pixel, not 4x4. At 4x4 the coverage of a pixel can only take the
     * seventeen values k/16, and the marks here are sub-pixel -- a cell at Wang's density
     * has a radius of 0.6 px, so most of the pixels it touches have true coverage well
     * under 1/16 and every one of them rounds up to it. Multiplied by 273,508 cells that
     * turned a faint rim into a visible one: ink coverage came out at 44% against
     * matplotlib's 35% on the isolate view, where nearly every cell is drawn in the
     * low-contrast context grey and a small alpha error decides whether a pixel reads as
     * drawn at all. The stamps are cached per (radius, sub-pixel bucket), so the finer
     * grid is paid for a few dozen times in a render, not 273,508 times. */
    const SS = 16, inv = 1 / SS, base = inv / 2;
    for (let py = 0; py < size; py++) {
      for (let px = 0; px < size; px++) {
        let hit = 0;
        for (let sy = 0; sy < SS; sy++) {
          for (let sx = 0; sx < SS; sx++) {
            const dx = (px - half) - 0.5 + base + sx * inv - fx;
            const dy = (py - half) - 0.5 + base + sy * inv - fy;
            if (dx * dx + dy * dy <= rr * rr) hit++;
          }
        }
        a[py * size + px] = hit / (SS * SS);
      }
    }
    const st = { half: half, size: size, a: a };
    stampCache.set(key, st);
    return st;
  }

  /* Which cached stamp a mark at `c` wants, and which pixel it lands in. */
  function bucket(c) {
    const i = Math.floor(c);
    let b = Math.floor((c - i) * SUBPIX);
    if (b < 0) b = 0; else if (b >= SUBPIX) b = SUBPIX - 1;
    return [i, b];
  }

  function blit(buf, W, H, ix, iy, st, rgb, alphaScale) {
    const x0 = ix - st.half, y0 = iy - st.half;
    const R = rgb[0], G = rgb[1], B = rgb[2];
    for (let py = 0; py < st.size; py++) {
      const yy = y0 + py;
      if (yy < 0 || yy >= H) continue;
      const rowA = py * st.size, rowB = yy * W;
      for (let px = 0; px < st.size; px++) {
        const a = st.a[rowA + px] * alphaScale;
        if (a <= 0) continue;
        const xx = x0 + px;
        if (xx < 0 || xx >= W) continue;
        const o = (rowB + xx) * 4;
        const ia = 1 - a;
        buf[o] = buf[o] * ia + R * a;
        buf[o + 1] = buf[o + 1] * ia + G * a;
        buf[o + 2] = buf[o + 2] * ia + B * a;
      }
    }
  }

  /* Deterministic shuffle, as in api.js -- see the note there on why a fixed permutation
   * is the requirement and numpy's particular one is not. */
  function shuffled(idx, seed) {
    let a = (seed >>> 0);
    const rnd = function () {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    for (let i = idx.length - 1; i > 0; i--) {
      const j = Math.floor(rnd() * (i + 1));
      const t = idx[i]; idx[i] = idx[j]; idx[j] = t;
    }
    return idx;
  }

  /* The drawing itself, into a plain RGBA buffer, with no DOM anywhere in it.
   *
   * Split out from CV.background so it can be run and checked off-browser. This is the
   * one module whose output is a picture rather than a number, and a picture is exactly
   * what could be silently wrong -- the y axis flipped, the extent off by a margin, the
   * context and selected sets the wrong way round -- while every other test passes. With
   * the pixels reachable from a script, they can be held against the matplotlib render
   * that this replaces (tests/compare_raster.js). `opts` matches the query string the
   * served /api/bg takes, so the caller builds one description and uses it as both the
   * cache key and the argument list. */
  CV.rasterBuffer = function (b, opts) {
    const mode = opts.mode || 'celltype';
    if (['celltype', 'profile', 'plain'].indexOf(mode) < 0) {
      throw new Error('unknown mode ' + mode);
    }
    const select = (opts.select === undefined || opts.select === null) ? null : +opts.select;
    const profile = (opts.profile === undefined || opts.profile === null) ? null : +opts.profile;
    const mana = !!opts.mana;
    const width = Math.max(200, Math.min(2400, Math.round(opts.width || 1400)));

    const ext = b.meta.extent.cell;
    const [ex0, ex1, ey0, ey1] = ext;
    const aspect = (ey1 - ey0) / (ex1 - ex0);
    const height = Math.max(1, Math.round(width * aspect));

    const buf = new Float32Array(width * height * 4);

    const ui = b.colours.ui;
    const surf = hexToRgb(ui.surface);
    for (let i = 0; i < width * height; i++) {
      buf[i * 4] = surf[0]; buf[i * 4 + 1] = surf[1];
      buf[i * 4 + 2] = surf[2]; buf[i * 4 + 3] = 255;
    }

    // The axes fill the figure exactly (render.py's add_axes([0,0,1,1])) and the limits
    // are the extent, so data maps to pixels linearly with y flipped.
    const sx = width / (ex1 - ex0), sy = height / (ey1 - ey0);
    const PX = (x) => (x - ex0) * sx;
    const PY = (y) => height - (y - ey0) * sy;

    const n = b.n_cells;
    const s = pointSize(n, width);
    const ctxRgb = hexToRgb(ui.context);

    /* One mark, at its true sub-pixel position, in the given colour. `marker` picks
     * between matplotlib's '.' and 'o' sizing. */
    function mark(r, x, y, rgb, marker) {
      const px = PX(x), py = PY(y);
      const bx = bucket(px), by = bucket(py);
      blit(buf, width, height, bx[0], by[0],
           stamp(r, bx[1], by[1]), rgb, 1);
    }
    const rCtx = radiusPx(s, '.');

    let values = null, colourFor = null, classes = [];
    if (mode === 'profile') {
      values = b.c_profile;
      colourFor = (k) => b.profile_colours[k];
      classes = Object.keys(b.profile_colours).map(Number).sort((p, r) => p - r);
    } else if (mode === 'celltype') {
      values = b.c_type;
      colourFor = (k) => b.cell_type_colours[k];
      classes = [];
      for (let i = 0; i < b.cell_type_colours.length; i++) classes.push(i);
    }

    if (values === null) {
      for (let i = 0; i < n; i++) {
        mark(rCtx, b.cxy[2 * i], b.cxy[2 * i + 1], ctxRgb);
      }
    } else {
      /* ONE MASK. `ctx` decides what is context, and the coloured set is exactly its
       * complement -- not recomputed from `values`. Deriving the two independently is how
       * the profile filter was silently dropped once in the served version: the context
       * mask had the filter in it and the selection did not, so asking for one profile
       * returned the whole map coloured. */
      const isCtx = new Uint8Array(n);
      for (let i = 0; i < n; i++) {
        let c;
        if (select === null) c = (mode === 'celltype') ? (values[i] < 0) : (values[i] <= 0);
        else c = (values[i] !== select);
        if (profile !== null && b.c_profile[i] !== profile) c = true;
        isCtx[i] = c ? 1 : 0;
      }
      for (let i = 0; i < n; i++) {
        if (!isCtx[i]) continue;
        mark(rCtx, b.cxy[2 * i], b.cxy[2 * i + 1], ctxRgb);
      }

      const idx = [];
      for (let i = 0; i < n; i++) if (!isCtx[i]) idx.push(i);
      if (idx.length) {
        shuffled(idx, 0);
        // Under a filter the selection is figure and the rest is ground, so the selected
        // cells are drawn larger -- without it a diffuse profile lays a thin scatter over
        // 273k context points and is close to invisible.
        const selS = (select !== null || profile !== null) ? s * 2.6 : s;
        const rSel = radiusPx(selS, '.');
        const rgbOf = {};
        for (const k of classes) {
          const h = colourFor(k);
          if (h) rgbOf[k] = hexToRgb(h);
        }
        for (const i of idx) {
          const rgb = rgbOf[values[i]];
          if (!rgb) continue;
          mark(rSel, b.cxy[2 * i], b.cxy[2 * i + 1], rgb);
        }
      }
    }

    /* MANAscore last, so it is never buried. Larger, red, and with a white halo rather
     * than relying on hue alone: the categorical palettes already contain reds, so on a
     * dense map the colour by itself would not reliably separate a called cell. */
    if (mana && b.has_manascore) {
      // marker='o' here, not '.', so no point-size reduction -- and the white edge is a
      // stroke of `s * 0.16` POINTS, which straddles the marker's rim, so it extends half
      // its width beyond it.
      const manaS = s * 3.4;
      const rDot = radiusPx(manaS, 'o');
      /* THE STROKE IS CENTRED ON THE MARKER'S EDGE, so it covers the outer half of what
       * would otherwise be red: matplotlib fills the path to rDot and then strokes a
       * white line of `s * 0.16` points along that same boundary, leaving red only inside
       * rDot - lw/2. Painting red out to rDot -- which is what this did first -- makes the
       * red disc 46% larger in area than the one being replaced, and since the halo is
       * drawn under it the marker simply looks bigger and redder rather than wrong. */
      const lwPx = (s * 0.16) * PT_PX;
      const rOut = rDot + lwPx * 0.5;
      const rRed = Math.max(0.35, rDot - lwPx * 0.5);
      const white = [255, 255, 255];
      const red = hexToRgb(ui.mana);
      for (let i = 0; i < n; i++) {
        if (b.c_mana[i] !== 1) continue;
        mark(rOut, b.cxy[2 * i], b.cxy[2 * i + 1], white);
        mark(rRed, b.cxy[2 * i], b.cxy[2 * i + 1], red);
      }
    }

    return { width: width, height: height, buf: buf };
  };

  /* The browser-facing wrapper: the same pixels, in a canvas. */
  CV.background = function (b, opts) {
    const r = CV.rasterBuffer(b, opts);
    const cv = document.createElement('canvas');
    cv.width = r.width; cv.height = r.height;
    const g = cv.getContext('2d');
    const img = g.createImageData(r.width, r.height);
    const out = img.data;
    for (let i = 0; i < r.width * r.height * 4; i += 4) {
      out[i] = r.buf[i]; out[i + 1] = r.buf[i + 1];
      out[i + 2] = r.buf[i + 2]; out[i + 3] = 255;
    }
    g.putImageData(img, 0, 0);
    return cv;
  };

  CV._raster = { pointSize, radiusPx, stamp, hexToRgb, shuffled };
}(window.CV));
