
#!/bin/bash
set -euo pipefail

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate
outdir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/vcf_chroms_MNP_annotate_samples

curr_pgm=$0
curr_pgm_base=$(basename $curr_pgm)
curr_pgm_base="${curr_pgm_base%.sh}"
out_manifest=./manifests/manifest."${curr_pgm_base}".txt

:>"${out_manifest}"

for infile in "${indir}"/"${batch}".gatk_vqsr_and_platypus_mnp.*.annovar_clinvarDATE.vep_reformatted.extras.vcf; do

  infile_basename=$(basename $infile)
  IFS='.' read -r -a array <<< "$infile_basename"
  prefix="${array[0]}"
  suffix="${infile_basename/$prefix/}"

  outfile_basename="ExtractSamples${suffix}"

  echo -e "${infile}\t${outdir}\t${outfile_basename}" >> "${out_manifest}"

done

echo 'out_manifest:' $out_manifest

./extract_multiple_samples_from_jointcalled_vcf_by_chrom_SubmitJobs.sh $out_manifest

