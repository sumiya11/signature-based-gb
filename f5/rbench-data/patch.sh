#!/usr/bin/env bash

awk '
BEGIN {all=1}
/share/ {all=0}
(all==1) {print; next}
/gb :=/ {print $0 "\n"}
/end/ {print}
' $1 > tmp && mv tmp $1
