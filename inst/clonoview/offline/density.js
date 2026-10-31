/* Density contours, ported from clonoview/density.py.
 *
 * Read that file first -- the statistics, the thresholds and the reasons for them live
 * there and are not repeated here. This is a translation, and it is deliberately a
 * literal one: the same bandwidth rule, the same bandwidth-sized grid, the same
 * highest-density-region levels, the same three refusals. Where a numpy or scipy call had
 * a behaviour that is easy to get subtly wrong (round-half-to-even, linear quantile
 * interpolation, the exact Gaussian kernel scipy builds) that behaviour is reimplemented
 * here rather than approximated, because the self-check compares the two and a silent
 * half-bin drift would be indistinguishable from a bug in the port.
 *
 * THE ONE PLACE THIS CANNOT MATCH PYTHON EXACTLY is the final step. matplotlib's contour
 * generator emits the same *geometry* as any marching-squares implementation, but the
 * order it walks cells in, where it starts a closed loop, and which way round it winds
 * are internal to it. Vertex-by-vertex comparison against it is therefore not a
 * meaningful test, and selftest.html compares the invariants instead: the level values,
 * how many polygons, their total enclosed area, their total perimeter, and their bounding
 * box. Two contour implementations that agree on all of those over hundreds of clones are
 * drawing the same rings.
 */

window.CV = window.CV || {};

(function (CV) {
  'use strict';

  const MIN_FOR_CONTOUR = 30;
  const LEVELS = [0.9, 0.5, 0.25];
  const GRID = 140;
  const TARGET_SIGMA = 8.0;

  /* numpy rounds halves to even; JavaScript's Math.round rounds halves up. This decides
   * the grid size, so on a tie the two would build grids of different sizes and produce
   * genuinely different contours -- a loud failure with an obscure cause. */
  function rint(v) {
    const f = Math.floor(v);
    const d = v - f;
    if (d > 0.5) return f + 1;
    if (d < 0.5) return f;
    return (f % 2 === 0) ? f : f + 1;
  }

  function stdDdof1(vals) {
    const n = vals.length;
    if (n < 2) return 0;
    let m = 0;
    for (let i = 0; i < n; i++) m += vals[i];
    m /= n;
    let s = 0;
    for (let i = 0; i < n; i++) { const d = vals[i] - m; s += d * d; }
    return Math.sqrt(s / (n - 1));
  }

  /* np.quantile's default 'linear' method: position q*(m-1) in the sorted values, with
   * linear interpolation between the two neighbours. Several other conventions exist and
   * they disagree by up to one whole order statistic on small samples, which is exactly
   * the regime the smaller clones sit in. */
  function quantileLinear(sorted, q) {
    const m = sorted.length;
    if (m === 0) return NaN;
    if (m === 1) return sorted[0];
    const pos = q * (m - 1);
    const lo = Math.floor(pos);
    const frac = pos - lo;
    if (lo >= m - 1) return sorted[m - 1];
    return sorted[lo] + frac * (sorted[lo + 1] - sorted[lo]);
  }

  /* scipy.ndimage.gaussian_filter's kernel, built the way scipy builds it: radius from
   * the default truncate of 4.0 as int(4*sigma + 0.5), an unnormalised Gaussian over
   * -radius..radius, then divided by its own sum. Normalising the *truncated* kernel
   * rather than using the analytic constant is why the result is not simply a sampled
   * Gaussian, and it matters most at the small sigmas the compact clones produce. */
  function gaussianKernel1d(sigma) {
    const radius = Math.floor(4.0 * sigma + 0.5);
    const k = new Float64Array(2 * radius + 1);
    let sum = 0;
    for (let i = -radius; i <= radius; i++) {
      const v = Math.exp(-0.5 * (i * i) / (sigma * sigma));
      k[i + radius] = v;
      sum += v;
    }
    for (let i = 0; i < k.length; i++) k[i] /= sum;
    return { k: k, radius: radius };
  }

  /* Separable correlation with mode='constant', cval=0: everything outside the array
   * counts as zero. For a density built from a histogram that is the right boundary --
   * there are no cells outside the grid, so there is no mass there. */
  function gaussianFilter2d(Z, ny, nx, sigma0, sigma1) {
    const g0 = gaussianKernel1d(sigma0);
    const g1 = gaussianKernel1d(sigma1);
    const tmp = new Float64Array(ny * nx);
    // axis 0 (down columns)
    for (let j = 0; j < nx; j++) {
      for (let i = 0; i < ny; i++) {
        let acc = 0;
        for (let t = -g0.radius; t <= g0.radius; t++) {
          const ii = i + t;
          if (ii < 0 || ii >= ny) continue;
          acc += Z[ii * nx + j] * g0.k[t + g0.radius];
        }
        tmp[i * nx + j] = acc;
      }
    }
    const out = new Float64Array(ny * nx);
    // axis 1 (along rows)
    for (let i = 0; i < ny; i++) {
      const row = i * nx;
      for (let j = 0; j < nx; j++) {
        let acc = 0;
        for (let t = -g1.radius; t <= g1.radius; t++) {
          const jj = j + t;
          if (jj < 0 || jj >= nx) continue;
          acc += tmp[row + jj] * g1.k[t + g1.radius];
        }
        out[row + j] = acc;
      }
    }
    return out;
  }

  /* Marching squares. Segments first, chained into polylines second.
   *
   * The ambiguous case -- two diagonally opposite corners above the level and two below,
   * which can be read as one saddle or two separate corners -- is resolved by the mean of
   * the four corners, which is what matplotlib does. Choosing the other reading would
   * join two rings that should be separate, or split one that should not be, so this is
   * not a detail that can be left to whichever branch happens to be written first.
   */
  function marchingSquares(Z, gx, gy, nx, ny, level) {
    const segs = [];
    const at = (i, j) => Z[i * nx + j];
    const interp = (v0, v1, c0, c1) => {
      const d = v1 - v0;
      if (d === 0) return c0;
      return c0 + (level - v0) / d * (c1 - c0);
    };
    for (let i = 0; i < ny - 1; i++) {
      for (let j = 0; j < nx - 1; j++) {
        const z00 = at(i, j), z10 = at(i, j + 1);
        const z11 = at(i + 1, j + 1), z01 = at(i + 1, j);
        let idx = 0;
        if (z00 > level) idx |= 1;
        if (z10 > level) idx |= 2;
        if (z11 > level) idx |= 4;
        if (z01 > level) idx |= 8;
        if (idx === 0 || idx === 15) continue;
        const x0 = gx[j], x1 = gx[j + 1], y0 = gy[i], y1 = gy[i + 1];
        // Edge crossing points: bottom, right, top, left of the cell.
        const B = [interp(z00, z10, x0, x1), y0];
        const R = [x1, interp(z10, z11, y0, y1)];
        const T = [interp(z01, z11, x0, x1), y1];
        const L = [x0, interp(z00, z01, y0, y1)];
        const push = (a, b) => segs.push([a[0], a[1], b[0], b[1]]);
        switch (idx) {
          case 1: case 14: push(L, B); break;
          case 2: case 13: push(B, R); break;
          case 3: case 12: push(L, R); break;
          case 4: case 11: push(R, T); break;
          case 6: case 9: push(B, T); break;
          case 7: case 8: push(L, T); break;
          case 5: case 10: {
            const mid = (z00 + z10 + z11 + z01) / 4;
            const centreAbove = mid > level;
            if (idx === 5) {                     // corners 00 and 11 above
              if (centreAbove) { push(L, T); push(B, R); }
              else { push(L, B); push(R, T); }
            } else {                             // corners 10 and 01 above
              if (centreAbove) { push(L, B); push(R, T); }
              else { push(L, T); push(B, R); }
            }
            break;
          }
          default: break;
        }
      }
    }
    return chain(segs);
  }

  /* Join segments end to end into polylines, closing loops where they close.
   *
   * Endpoints are matched through a key quantised to a grid far finer than the data, not
   * by exact float equality: the same crossing point computed from the two cells that
   * share an edge is arithmetically identical here (both evaluate the same interp on the
   * same two corner values) but relying on that would make the chaining hostage to a
   * later refactor that reorders the arithmetic. Unmatched segments are still returned,
   * as one-segment paths, rather than dropped -- a ring that fails to close is a visible
   * gap the reader can see, whereas a dropped one silently shrinks the contour.
   */
  function chain(segs) {
    const KEY = 1e9;
    const key = (x, y) => Math.round(x * KEY) + ',' + Math.round(y * KEY);
    const ends = new Map();
    const used = new Uint8Array(segs.length);
    for (let i = 0; i < segs.length; i++) {
      const s = segs[i];
      for (const k of [key(s[0], s[1]), key(s[2], s[3])]) {
        if (!ends.has(k)) ends.set(k, []);
        ends.get(k).push(i);
      }
    }
    const paths = [];
    for (let i = 0; i < segs.length; i++) {
      if (used[i]) continue;
      used[i] = 1;
      const s = segs[i];
      const pts = [[s[0], s[1]], [s[2], s[3]]];
      // Walk forward from the tail, then backward from the head.
      for (let dir = 0; dir < 2; dir++) {
        for (;;) {
          const p = dir === 0 ? pts[pts.length - 1] : pts[0];
          const cand = ends.get(key(p[0], p[1]));
          if (!cand) break;
          let nxt = -1;
          for (const c of cand) if (!used[c]) { nxt = c; break; }
          if (nxt < 0) break;
          used[nxt] = 1;
          const t = segs[nxt];
          const head = (Math.round(t[0] * KEY) === Math.round(p[0] * KEY) &&
                        Math.round(t[1] * KEY) === Math.round(p[1] * KEY));
          const other = head ? [t[2], t[3]] : [t[0], t[1]];
          if (dir === 0) pts.push(other); else pts.unshift(other);
        }
      }
      paths.push(pts);
    }
    return paths;
  }

  /* xs, ys: parallel arrays of the cells' coordinates. extent: [x0,x1,y0,y1] of the
   * whole embedding, which caps the padded grid so a clone at the edge of the map does
   * not get a grid running off it. */
  CV.contours = function (xs, ys, extent, opts) {
    opts = opts || {};
    const levels = opts.levels || LEVELS;
    const grid = opts.grid || GRID;
    const minCells = opts.minCells === undefined ? MIN_FOR_CONTOUR : opts.minCells;
    const n = xs.length;
    const out = { n: n, levels: [], rings: [], reason: null };
    if (n < minCells) {
      out.reason = n + ' cells — below the ' + minCells + '-cell minimum for a ' +
                   'density estimate; cells shown individually';
      return out;
    }

    /* opts.sd exists for the comparison harness only. numpy computes the standard
     * deviation of a float32 array in float32, which is about 1e-8 less accurate
     * relatively than the float64 done here -- and on a small clone that tiny difference
     * is amplified, because the contour levels are read from the density AT the cells
     * through a grid lookup, so a coordinate shift of 1e-8 can move a cell across a bin
     * edge and change an order statistic outright. Being able to inject Python's own
     * standard deviation is what separates "the port disagrees" from "the two arithmetics
     * disagree", and the answer decides whether there is a bug to fix. */
    const sdx = opts.sd ? opts.sd[0] : stdDdof1(xs);
    const sdy = opts.sd ? opts.sd[1] : stdDdof1(ys);
    if (!isFinite(sdx) || !isFinite(sdy) || (sdx <= 0 && sdy <= 0)) {
      out.reason = 'cells occupy no area; shown individually';
      return out;
    }
    // np.where(sd > 0, sd, np.max(sd)): a clone lying on a perfectly horizontal or
    // vertical line has one zero standard deviation, and borrows the other's.
    const mx = Math.max(sdx, sdy);
    const f = Math.pow(n, -1.0 / 6.0);
    const bwx = (sdx > 0 ? sdx : mx) * f;
    const bwy = (sdy > 0 ? sdy : mx) * f;

    let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
    for (let i = 0; i < n; i++) {
      if (xs[i] < xmin) xmin = xs[i];
      if (xs[i] > xmax) xmax = xs[i];
      if (ys[i] < ymin) ymin = ys[i];
      if (ys[i] > ymax) ymax = ys[i];
    }
    const x0 = Math.max(extent[0], xmin - 3 * bwx);
    const x1 = Math.min(extent[1], xmax + 3 * bwx);
    const y0 = Math.max(extent[2], ymin - 3 * bwy);
    const y1 = Math.min(extent[3], ymax + 3 * bwy);
    if (!(x1 > x0 && y1 > y0)) {
      out.reason = 'cells occupy no area; shown individually';
      return out;
    }

    const nx = Math.min(grid, Math.max(48, rint(TARGET_SIGMA * (x1 - x0) / bwx)));
    const ny = Math.min(grid, Math.max(48, rint(TARGET_SIGMA * (y1 - y0) / bwy)));

    const gx = new Float64Array(nx), gy = new Float64Array(ny);
    for (let i = 0; i < nx; i++) gx[i] = x0 + (x1 - x0) * i / (nx - 1);
    for (let i = 0; i < ny; i++) gy[i] = y0 + (y1 - y0) * i / (ny - 1);

    /* np.histogram2d over nx bins across [x0,x1] and ny across [y0,y1]. The bin count is
     * the number of GRID POINTS, so the bin edges and the grid coordinates are offset by
     * half a bin from each other. That is what density.py does and it is reproduced, not
     * corrected: changing it would move every existing ring by half a bin. */
    const H = new Float64Array(ny * nx);
    const dx = (x1 - x0) / nx, dy = (y1 - y0) / ny;
    for (let i = 0; i < n; i++) {
      const xv = xs[i], yv = ys[i];
      if (xv < x0 || xv > x1 || yv < y0 || yv > y1) continue;   // outside range: dropped
      let bx = Math.floor((xv - x0) / dx);
      let by = Math.floor((yv - y0) / dy);
      if (bx >= nx) bx = nx - 1;                                // last bin closed on the right
      if (by >= ny) by = ny - 1;
      if (bx < 0) bx = 0;
      if (by < 0) by = 0;
      H[by * nx + bx] += 1;
    }

    const sigma0 = Math.max(0.5, bwy / (gy[1] - gy[0]));
    const sigma1 = Math.max(0.5, bwx / (gx[1] - gx[0]));
    const Z = gaussianFilter2d(H, ny, nx, sigma0, sigma1);

    let zmax = -Infinity, zmin = Infinity;
    for (let i = 0; i < Z.length; i++) {
      if (Z[i] > zmax) zmax = Z[i];
      if (Z[i] < zmin) zmin = Z[i];
    }
    if (zmax <= 0) {
      out.reason = 'density is flat at this resolution; cells shown individually';
      return out;
    }

    /* The level for fraction f is read off the density AT the cells, on the same grid
     * the contours come from, so the ring and the number it is labelled with cannot
     * disagree. np.searchsorted(side='left') minus one, clipped. */
    const atPoints = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let jx = searchsortedLeft(gx, xs[i]) - 1;
      let jy = searchsortedLeft(gy, ys[i]) - 1;
      if (jx < 0) jx = 0; else if (jx > nx - 1) jx = nx - 1;
      if (jy < 0) jy = 0; else if (jy > ny - 1) jy = ny - 1;
      atPoints[i] = Z[jy * nx + jx];
    }
    const sorted = Float64Array.from(atPoints).sort();

    const want = [];
    for (const fr of levels) {
      const lv = quantileLinear(sorted, 1.0 - fr);
      if (lv > zmin) want.push([fr, lv]);
    }
    if (!want.length) {
      out.reason = 'density is flat at this resolution; cells shown individually';
      return out;
    }
    want.sort((a, b) => a[1] - b[1]);
    const keep = [want[0]];
    for (let i = 1; i < want.length; i++) {
      if (want[i][1] > want[i - 1][1]) keep.push(want[i]);
    }

    /* An inspection hook, off unless something sets it. The contour is four stages deep
     * -- bin, smooth, pick levels, trace -- and a disagreement with Python at the end
     * says nothing about which stage produced it. Being able to lift the intermediate
     * grid out is what turned a 0.8% area difference from a mystery into a located
     * cause; it costs one branch and it is the difference between checking this port and
     * arguing about it. */
    if (CV.INSPECT) CV.INSPECT({ nx, ny, gx, gy, Z, keep, sorted, zmin, zmax });

    for (const [fr, lv] of keep) {
      const paths = marchingSquares(Z, gx, gy, nx, ny, lv);
      const rings = [];
      for (const p of paths) {
        if (p.length < 3) continue;      // matplotlib's own >= 3 vertices rule
        const flat = [];
        for (const pt of p) { flat.push(pt[0], pt[1]); }
        rings.push(flat);
      }
      out.levels.push(Math.round(fr * 1000) / 1000);
      out.rings.push(rings);
    }
    return out;
  };

  function searchsortedLeft(arr, v) {
    let lo = 0, hi = arr.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (arr[mid] < v) lo = mid + 1; else hi = mid;
    }
    return lo;
  }

  CV.MIN_FOR_CONTOUR = MIN_FOR_CONTOUR;
  CV._density = { rint, quantileLinear, gaussianKernel1d, gaussianFilter2d,
                  marchingSquares, searchsortedLeft };
}(window.CV));
