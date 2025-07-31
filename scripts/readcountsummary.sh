#!/bin/bash

for fn in $*;
do
   printf "$fn\t"
   cat "$fn"
done