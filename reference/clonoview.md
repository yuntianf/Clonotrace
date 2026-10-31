# clonoview

Serve a prepared dataset in the ClonoView viewer.

## Usage

``` r
clonoview(dataset = NULL, port = NULL, open = interactive(), python = NULL)
```

## Arguments

- dataset:

  Dataset to open on load. With more than one prepared and none named,
  the viewer opens on a chooser.

- port:

  Port to serve on. Defaults to the first free port from 8770.

- open:

  Whether to open a browser. Defaults to TRUE in an interactive session
  on a machine that has one.

- python:

  Python interpreter to use. Defaults to \[clonoview_python()\].

## Value

The URL, invisibly.

## Details

Starts the viewer in the background and returns, so the R session stays
usable while the browser is open. Stop it with \[clonoview_stop()\].

The server binds to the loopback address, which means it is reachable
from the machine it runs on and from nowhere else. On a laptop that is
all that is needed and the browser opens straight onto it. On a cluster
the way in is an SSH tunnel, and the command to open one is printed;
binding to the node's own address instead would put unpublished data on
every host that can route to it, for no gain over the tunnel.

## See also

\[clonoview_bundle()\], \[clonoview_stop()\]
