#!/bin/bash

cat $1 | grep ">" | perl -lape 's/>(\S+).*$/$1/g' | sort
