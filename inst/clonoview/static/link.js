/* The link between the viewer tab and the table tab.
 *
 * WHY THIS IS NOT JUST A BroadcastChannel ANY MORE. A BroadcastChannel is scoped to an
 * origin, and a page opened from file:// -- which is how the offline export runs -- is
 * given an OPAQUE origin that matches nothing, not even another copy of itself in the
 * next tab. So in the export the two tabs would each open a channel, each post into it,
 * and neither would ever hear the other: the tables would still render, clicking a row
 * would still look like it did something, and the viewer would simply never change. That
 * is the worst failure shape there is, and it is invisible until someone clicks.
 *
 * window.postMessage is the exception. It is designed to cross origins -- that is its
 * entire purpose -- so it works between an opener and the window it opened even when both
 * are opaque file:// pages. It needs a handle rather than a name: the viewer keeps what
 * window.open() returned, and the table uses window.opener.
 *
 * BOTH ARE USED, AND DUPLICATES ARE DROPPED. postMessage needs the opener relationship,
 * which is missing if the table tab is reloaded from a bookmark or opened by hand;
 * BroadcastChannel does not, but only works when served. Rather than decide which
 * situation this is -- and be wrong in the one case nobody tested -- every message goes
 * out on both and arrives tagged, and the receiver ignores an id it has already seen.
 * Sending twice is cheap; a silent link is not.
 */

window.LINK = (function () {
  'use strict';

  const NAME = 'clonoview';
  const TAG = Math.random().toString(36).slice(2) + '-' + (window.opener ? 'tab' : 'main');
  let seq = 0;

  const chan = ('BroadcastChannel' in window) ? new BroadcastChannel(NAME) : null;
  let child = null;                      // set by the viewer when it opens the table tab
  const handlers = [];
  const seen = [];                       // recent message ids, newest last
  const SEEN_MAX = 64;

  function peers() {
    const out = [];
    if (child && !child.closed) out.push(child);
    if (window.opener && !window.opener.closed) out.push(window.opener);
    return out;
  }

  function send(msg) {
    const wrapped = Object.assign({}, msg, { _id: TAG + ':' + (++seq), _from: TAG });
    // targetOrigin '*': the receiver is a file:// page with an opaque origin, so there is
    // no origin string that would match it. Nothing secret travels over this link -- it
    // carries a lineage index and a list of patient indices, both of which the receiving
    // page already has in full -- and a page can only be a peer here by having been
    // opened by this one or having opened it.
    for (const w of peers()) {
      try { w.postMessage(wrapped, '*'); } catch (e) { /* window went away */ }
    }
    if (chan) {
      try { chan.postMessage(wrapped); } catch (e) { /* channel closed */ }
    }
  }

  function deliver(m) {
    if (!m || typeof m !== 'object' || !m._id) return;
    if (m._from === TAG) return;                 // our own, echoed back
    if (seen.indexOf(m._id) >= 0) return;        // already arrived on the other transport
    seen.push(m._id);
    if (seen.length > SEEN_MAX) seen.shift();
    for (const h of handlers) {
      try { h(m); } catch (e) { /* one handler must not stop the others */ }
    }
  }

  window.addEventListener('message', (ev) => deliver(ev.data));
  if (chan) chan.onmessage = (ev) => deliver(ev.data);

  return {
    send: send,
    onMessage: function (fn) { handlers.push(fn); },
    attachChild: function (w) { child = w; },
    /* Whether a message posted now could actually reach anyone. The table tab checks this
     * to decide what to tell the reader, and it has to be honest: "linked" when nothing
     * is listening is exactly the lie this file exists to prevent. */
    connected: function () { return peers().length > 0 || !!chan; },
    describe: function () {
      const p = peers().length;
      if (p && chan) return 'linked to the viewer tab';
      if (p) return 'linked to the viewer tab';
      if (chan) return 'linked to the viewer tab';
      return 'NOT linked — this browser cannot connect the two tabs';
    },
  };
}());
