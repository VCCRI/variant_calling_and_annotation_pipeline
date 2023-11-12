#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/manta_annotate
database_file=/my/directory/variant_calling_and_annotation_pipeline/code/databases/database.manta.allVars.txt

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/*.manta.CDSx.Gene.bndX.gCDSx.gGene.gBndX.EHRF.gSV.DGV.Dup.dbsc.TSS.pext.Idel.Ishort.gIshort.hot.collapsed.tsv; do

  sample=$(basename $infile)
  IFS='.' read -r -a array <<< "$sample"
  sample="${array[0]}"

  exclude_samples="${sample}"
  chr1to3=${sample:0:3}
  if [[ $chr1to3 == "FAM" ]]; then
    IFS='_' read -r -a array <<< "$sample"
    exclude_samples=${sample//\_/\|}
    #sample="${array[1]}"
  fi

  infile_basename=$(basename $infile)
  outfile_basename="${infile_basename%.txt}"
  outfile_basename="${outfile_basename%.tsv}"
  outfile_basename="${outfile_basename%.collapsed}"
  outfile_basename="${outfile_basename}".compareSamples.tsv

  include_low_qual_flag='.'

  echo -e "${sample}\t${infile}\t${outdir}\t${outfile_basename}\t${database_file}\t${include_low_qual_flag}\t${exclude_samples}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./manta_3b_compare_strucvars_with_other_samples_SubmitJobs.sh $out_manifest

