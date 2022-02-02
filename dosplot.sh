#!/bin/bash

scp wkdos/format.dat bern:~/as/wkdos/format.dat && \
ssh bern 'ssh bern6 "cd as/wkdos && ./dosout"' && \
scp bern:~/as/wkdos/dos.ai dos.ai && \
open dos.ai
