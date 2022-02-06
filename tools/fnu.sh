#!/bin/sh
find . -type f -name "fl*.dat" | xargs grep -in "umake"
