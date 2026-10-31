# clonoview_stop

Stop the running ClonoView server.

## Usage

``` r
clonoview_stop()
```

## Value

TRUE if a server was stopped, FALSE if none was running, invisibly.

## Details

Signals one process, named by the pidfile the server wrote once it had
bound. It does not match on the command line: a pattern like \`pkill -f
clonoview\` also matches the process doing the matching, and process ids
are reused, so on a shared machine that can signal somebody else's job.
