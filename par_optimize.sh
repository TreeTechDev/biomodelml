#!/bin/bash


if [ -z "$1" ]
  then
    CORES=1
  else
    CORES=$1
fi


for i in {1..$CORES}
do
	make optimize &
done