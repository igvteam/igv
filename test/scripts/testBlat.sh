#!/bin/bash
wget -q -O - --no-check-certificate "https://genome.ucsc.edu/cgi-bin/hgBlat?userSeq=$1&type=DNA&db=$2&output=json"
