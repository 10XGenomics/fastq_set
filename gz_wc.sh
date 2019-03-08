#!/bin/bash

gzip -c -d $1 | wc -l
