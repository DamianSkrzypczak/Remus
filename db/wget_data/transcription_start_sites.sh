#!/usr/bin/env bash
mkdir transcription_start_sites_data
wget -O promoter_data.bed 'http://promoter.binf.ku.dk/viewer.php?match=and&sort-by=donotsort&end-site=249250621&start-site=1&chr-number=ALL&toggle=basic&return=download'
mv promoter_data.bed transcription_start_sites_data/