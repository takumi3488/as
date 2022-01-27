#!/bin/sh
scp bern:~/as/${@:$#:1} $@ && code ${@:$#:1}
