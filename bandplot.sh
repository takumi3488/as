#!/bin/bash

scp wkplot/fort.3 bern:~/as/wkplot/fort.3 && \
ssh bern 'ssh bern6 "cd as/wkplot && bandplot && pig2ps"' && \
scp bern:~/as/wkplot/fort.50 band.ai && \
open band.ai
