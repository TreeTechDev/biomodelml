#!/usr/bin/env bash
awk -v OFS='\t' {'print $1,$2'} $1 > $2