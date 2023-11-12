#!/bin/bash

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate/*.gridss.nofilter_svtypes.tsv

final_db=databases/database.gridss.allVars.include_low_qual.txt
old_db=databases/database.gridss.allVars.include_low_qual.OLD.txt

if [[ -f "$final_db" ]]; then
  mv $final_db $old_db
fi

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
tmp_db=databases/temp."${curr_pgm_base}".txt

:>$tmp_db
for infile in $indir; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  sample="${array[0]}"
  IFS='_' read -r -a array <<< "$sample"
  sample="${array[0]}"

  echo -e "${sample}\t${infile}" >> $tmp_db

done

if [[ -f "$old_db" ]]; then
  cat $tmp_db $old_db | sort | uniq > $final_db
else
  cat $tmp_db | sort | uniq > $final_db
fi

rm $tmp_db
