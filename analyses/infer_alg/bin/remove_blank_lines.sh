#!/usr/bin/env bash

file="$1"

awk -i inplace '$2 != ""' "$file"