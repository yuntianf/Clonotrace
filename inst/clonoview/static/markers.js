/* Marker shapes: the second categorical channel, shared by the viewer and the tables.
 *
 * WHY SHAPE AND NOT COLOUR. Fill colour already means the clone unit's profile, in
 * every panel of the viewer. Telling several hand-picked clones apart by colour would
 * therefore break the reading of the thing the viewer exists to show. Shape is
 * orthogonal: two markers of the same colour still agree about the profile and
 * disagree about which clone they belong to.
 *
 * WHY THIS IS ITS OWN FILE. The clone map draws these on a canvas (app.js) and the
 * clone table draws them as SVG beside a row (table.js), on two separate pages. Two
 * copies of the geometry would be two definitions of "which shape is clone 3", and
 * they would drift. There is one definition here and both pages read it.
 */
window.MARKERS = (function () {
  'use strict';

  const SHAPES = ['circle', 'square', 'triangle', 'diamond',
                  'triangleDown', 'plus', 'triangleLeft', 'triangleRight'];

  /* Beyond this many hand-picked clones, shape stops being readable, and repeating a
   * shape would mean two different clones look identical. Callers stop marking instead
   * and say so -- an ambiguous key is worse than an honest absence. */
  const MAX = SHAPES.length;

  /* Builds the outline on any object with the canvas path methods. Each shape is drawn
   * to roughly the area of the circle of the same r, so no clone looks bigger or more
   * important than another just by being a square. */
  function path(g, shape, px, py, r) {
    g.beginPath();
    switch (shape) {
      case 'square': {
        const a = r * 0.886;                    // equal area to the circle
        g.rect(px - a, py - a, 2 * a, 2 * a);
        break;
      }
      case 'diamond': {
        const a = r * 1.253;
        g.moveTo(px, py - a); g.lineTo(px + a, py);
        g.lineTo(px, py + a); g.lineTo(px - a, py);
        g.closePath();
        break;
      }
      case 'triangle': {
        const a = r * 1.5551;        // equilateral, equal area to the circle
        g.moveTo(px, py - a);
        g.lineTo(px + a * 0.866, py + a * 0.5);
        g.lineTo(px - a * 0.866, py + a * 0.5);
        g.closePath();
        break;
      }
      case 'triangleDown': {
        const a = r * 1.5551;        // equilateral, equal area to the circle
        g.moveTo(px, py + a);
        g.lineTo(px + a * 0.866, py - a * 0.5);
        g.lineTo(px - a * 0.866, py - a * 0.5);
        g.closePath();
        break;
      }
      case 'triangleLeft': {
        const a = r * 1.5551;        // equilateral, equal area to the circle
        g.moveTo(px - a, py);
        g.lineTo(px + a * 0.5, py - a * 0.866);
        g.lineTo(px + a * 0.5, py + a * 0.866);
        g.closePath();
        break;
      }
      case 'triangleRight': {
        const a = r * 1.5551;        // equilateral, equal area to the circle
        g.moveTo(px + a, py);
        g.lineTo(px - a * 0.5, py - a * 0.866);
        g.lineTo(px - a * 0.5, py + a * 0.866);
        g.closePath();
        break;
      }
      case 'plus': {
        const a = r * 1.0879, w = a * 0.42;   // equal area to the circle
        g.moveTo(px - w, py - a); g.lineTo(px + w, py - a);
        g.lineTo(px + w, py - w); g.lineTo(px + a, py - w);
        g.lineTo(px + a, py + w); g.lineTo(px + w, py + w);
        g.lineTo(px + w, py + a); g.lineTo(px - w, py + a);
        g.lineTo(px - w, py + w); g.lineTo(px - a, py + w);
        g.lineTo(px - a, py - w); g.lineTo(px - w, py - w);
        g.closePath();
        break;
      }
      default:
        g.arc(px, py, r, 0, 7);
    }
  }

  /* An <svg> glyph of one shape, for a legend or a table cell. It runs `path` against a
   * recording shim rather than reimplementing the outlines, so the glyph in the key
   * cannot drift away from the marker in the plot it explains. */
  function glyph(shape, size) {
    const s = size || 11, c = s / 2, r = s * 0.30;
    const pts = [];
    let round = false;
    const shim = {
      beginPath() {}, closePath() {},
      arc() { round = true; },
      moveTo(x, y) { pts.push([x, y]); },
      lineTo(x, y) { pts.push([x, y]); },
      rect(x, y, w, h) { pts.push([x, y], [x + w, y], [x + w, y + h], [x, y + h]); },
    };
    path(shim, shape, c, c, r);
    const open = `<svg class="mk" width="${s}" height="${s}" ` +
                 `viewBox="0 0 ${s} ${s}" aria-hidden="true">`;
    if (round || !pts.length) return open + `<circle cx="${c}" cy="${c}" r="${r}"/></svg>`;
    const d = pts.map(p => p[0].toFixed(2) + ',' + p[1].toFixed(2)).join(' ');
    return open + `<polygon points="${d}"/></svg>`;
  }

  /* Which shape the n-th hand-picked clone gets, or null when the selection is too
   * small to need marking or too large to mark honestly. */
  function shapeAt(i, count) {
    if (count < 2 || count > MAX) return null;
    return (i < 0 || i >= MAX) ? null : SHAPES[i];
  }

  return { SHAPES: SHAPES, MAX: MAX, path: path, glyph: glyph, shapeAt: shapeAt };
}());
