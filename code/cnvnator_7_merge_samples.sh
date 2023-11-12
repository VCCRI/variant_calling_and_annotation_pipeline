#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt # one line per family, or sample if sample is not part of a family
in_ped_file="${batch}".pedigree_file.ped # one line per samples, thus multiple lines per family

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/cnvnator_annotate

           infile_basename_suffix_CHD_singles=.cnvnator.CHD_955_genes.vpot_final_output_file.vpot_gt_0.txt
outfile_basename_suffix_CHD_singles=_MY_BATCH__cnvnator.CHD_955_genes.vpot_final_output_file.vpot_gt_0.txt

           infile_basename_suffix_DCM_singles=.cnvnator.DCM_909_genes.vpot_final_output_file.vpot_gt_0.txt
outfile_basename_suffix_DCM_singles=_MY_BATCH__cnvnator.DCM_909_genes.vpot_final_output_file.vpot_gt_0.txt

tmpdir="${PBS_JOBFS}"/tmp
mkdir -p "${tmpdir}"

cohorts=(CHD DCM)

# Call vpot for singletons, and also for famiies.
# Merged singletons report will get vpot-prioritised variants with vpot score > 0.
# Merged families report will get all vpot-prioritised variants because we will do a manual inheritance modelling.

for cohort in "${cohorts[@]}"; do

  echo 'Processing:' $cohort

  for singles_or_family in "singles"; do # cnvnator does not do joint calls, so there are no merged family results like there is for gridss and manta

    echo '   Sub-processing:' $singles_or_family

    infile_basename_suffix_varname=infile_basename_suffix_"${cohort}"_"${singles_or_family}"
    outfile_basename_suffix_varname=outfile_basename_suffix_"${cohort}"_"${singles_or_family}"

    infile_basename_suffix="${!infile_basename_suffix_varname}"
    outfile_basename_suffix="${!outfile_basename_suffix_varname}"

    outfile="${indir}"/"${cohort}_${singles_or_family}_${outfile_basename_suffix}"
    string_of_files=''

    if [[ $infile_basename_suffix != "" ]]; then

      # read through the manifest list of families or singles, each one will have a structural variant file,
      # then see if it is in the cohort and singles_or_family that we are currently processing to merge into one file
      one_file=''
      while IFS= read -r inline; do

        IFS=$'\t' read -r -a array <<< "$inline"
        this_flagship="${array[1]}"
        this_proband="${array[3]}"
        this_cohort="${array[0]}"
        this_family="${array[2]}"
        if [[ $this_family == "." ]]; then
          # this manifest entry is for a singleton
          this_singles_or_family="singles"
        else
          # this manifest entry is for a family
          this_singles_or_family="family"
        fi

        if [[ "$this_cohort" == "$cohort" ]]; then

          num_array_elements="${#array[@]}"
          for (( i=3; i<$num_array_elements; i++ )); do

            # This family/sample is in this cohort. All samples, even those in a trio or family, are processed as singletons by cnvnator

            # get the structural variants file for this singleton

            this_sample="${array[i]}"
            infile=$(ls -1 "${indir}"/"${this_sample}""${infile_basename_suffix}" | head -n 1)
            if [[ -f "$infile" ]]; then
              echo 'Found singleton file for' $this_sample ':' $infile
            else
              echo "ERROR: Cannot find file for" $this_sample
            fi

            echo ''
            echo 'Input file that will be used:' $infile
            echo ''

            new_string_of_files="${string_of_files} ${infile}"
            string_of_files=$new_string_of_files
            one_file=$infile

          done
        fi
      done < "$in_manifest"

      if [[ $one_file != '' ]]; then # some files were found in this cohort and singles_or_family, so output a merged file of those files

        tmphdr="${tmpdir}"/cnvnator_7_merge_samples.temp_cnvnator_hdr.txt
        tmphdr2="${tmpdir}"/cnvnator_7_merge_samples.temp_cnvnator_hdr2.txt

        head -n 1 $one_file > $tmphdr
        last_col=$(sed -e 's/\t/\n/g' "$tmphdr" | tail -n 1)
        sed_string='s/'"${last_col}"'$/sample_field/'
        sed $sed_string $tmphdr > $tmphdr2
        echo 'cat' $string_of_files '| grep -P -v' '^Rank|^chrom' '| sort -k1nr -k2nr | cat' $tmphdr2 '- >' $outfile '|| true'
        cat $string_of_files | grep -P -v '^Rank|^chrom' | sort -k1nr -k2nr | cat $tmphdr2 - > $outfile || true

        echo ''
        echo '##################################################'
        echo 'outfile:' $outfile
        echo '##################################################'
        echo ''
      fi

    fi
  done
done

