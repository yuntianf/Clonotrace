/* Loading and decoding the packed datasets, in the browser, with no server.
 *
 * WHY A SCRIPT TAG AND NOT fetch(). A page opened from file:// is given an opaque origin,
 * and fetch/XMLHttpRequest are refused across one -- there is no flag, header or CORS
 * dance that changes it, because there is no server to send a header. A classic
 * <script src> IS permitted, so that is the transport: pack.py writes each dataset as a
 * .js file that assigns itself into window.CV_DATA, and loadDataset() appends a tag and
 * waits for onload. ES modules are NOT an option either -- module scripts are fetched
 * under CORS rules and fail from file:// exactly as fetch does.
 *
 * WHY IT IS LAZY. The four datasets are 14.7 MB of JavaScript between them. The picker
 * has to be populated before the user chooses, so the catalogue (index.js, a few KB)
 * loads at startup and a dataset's arrays load only when it is opened. Each is decoded
 * once and kept -- switching back is then instant, and 14.7 MB is small enough that
 * holding all four is not worth managing.
 */

window.CV = window.CV || {};

(function (CV) {
  'use strict';

  /* Typed arrays read the platform's byte order and offer no way to ask for the other
   * one, so a big-endian machine would decode every coordinate as garbage rather than
   * failing. pack.py writes little-endian explicitly; this checks the reader agrees.
   * Every platform a browser runs on today is little-endian, which is exactly why this
   * would be missed if it ever stopped being true. */
  const LITTLE_ENDIAN = (function () {
    const b = new ArrayBuffer(2);
    new Uint16Array(b)[0] = 1;
    return new Uint8Array(b)[0] === 1;
  })();

  const VIEWS = {
    f4: Float32Array, i4: Int32Array, i2: Int16Array,
    i1: Int8Array, u1: Uint8Array,
  };

  function decodeArray(spec) {
    if (!LITTLE_ENDIAN) {
      throw new Error('This machine is big-endian; the packed data is little-endian.');
    }
    const bin = atob(spec.b64);
    const n = bin.length;
    const bytes = new Uint8Array(n);
    // charCodeAt in a plain loop rather than Uint8Array.from(bin, c => c.charCodeAt(0)):
    // the callback form allocates an intermediate for each of several million characters
    // and is several times slower on the 5.7 MB dataset, which is felt as a pause on the
    // dataset switch.
    for (let i = 0; i < n; i++) bytes[i] = bin.charCodeAt(i);
    return new VIEWS[spec.t](bytes.buffer, 0, spec.n * spec.cols);
  }

  /* The decoded counterpart of bundle.py's Bundle. Same field names on purpose: the
   * route implementations in api.js are read side by side with server.py's Api methods,
   * and renaming things across that boundary is how the two quietly diverge. */
  function Bundle(raw) {
    this.name = raw.name;
    this.meta = raw.meta;
    this.names = raw.names;
    this.colours = raw.colours;
    this.has_manascore = raw.has_manascore;
    const a = raw.arrays;
    for (const k in a) this[k] = decodeArray(a[k]);
    this.n_cells = raw.n_cells;
    this.n_units = raw.n_units;
    this.visits = this.meta.visits.map((v) => v.code);
    this.profile_ids = this.meta.profiles.map((p) => p.id);
    this.cell_types = this.meta.cell_types.map((t) => t.name);
    // profile_colours is keyed by the profile's own id, which is not its index: a run
    // whose Leiden step produced profiles 1,2,3,6,7 has no profile 4 or 5, and indexing
    // the palette positionally would shift every colour after the gap.
    this.profile_colours = {};
    for (const k in raw.colours.profile) this.profile_colours[+k] = raw.colours.profile[k];
    this.cell_type_colours = [];
    for (const k in raw.colours.cell_type) this.cell_type_colours[+k] = raw.colours.cell_type[k];
  }

  Bundle.prototype.cellsOfLineage = function (li) {
    return [this.lin_off[li], this.lin_off[li + 1]];
  };

  const loaded = {};
  const pending = {};

  CV.loadDataset = function (name) {
    if (loaded[name]) return Promise.resolve(loaded[name]);
    if (pending[name]) return pending[name];
    pending[name] = new Promise(function (resolve, reject) {
      if (window.CV_DATA && window.CV_DATA[name]) {
        resolve(loaded[name] = new Bundle(window.CV_DATA[name]));
        return;
      }
      const s = document.createElement('script');
      s.src = 'data/' + name + '.js';
      s.onload = function () {
        if (!window.CV_DATA || !window.CV_DATA[name]) {
          reject(new Error('data/' + name + '.js loaded but defined nothing'));
          return;
        }
        try {
          resolve(loaded[name] = new Bundle(window.CV_DATA[name]));
        } catch (e) { reject(e); }
      };
      // A script tag's onerror carries no detail -- for a missing file, a directory that
      // did not survive the unzip, or a syntax error the browser reports the same empty
      // event. The message says what to check rather than what happened, because the
      // person reading it unzipped a file and double-clicked something.
      s.onerror = function () {
        reject(new Error('Could not load data/' + name + '.js. If ClonoView was ' +
                         'unzipped, check that the "data" folder sits next to ' +
                         'ClonoView.html and still contains ' + name + '.js.'));
      };
      document.head.appendChild(s);
    });
    return pending[name];
  };

  CV.loadedDataset = function (name) { return loaded[name] || null; };
  CV.Bundle = Bundle;
  CV.decodeArray = decodeArray;
}(window.CV));
