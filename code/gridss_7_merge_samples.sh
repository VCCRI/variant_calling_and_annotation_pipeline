#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"

in_manifest="${batch}".one_line_per_family.txt # one line per family, or sample if sample is not part of a family
in_ped_file="${batch}".pedigree_file.ped # one line per samples, thus multiple lines per family

indir=/my/directory/variant_calling_and_annotation_pipeline/working_directory/gridss_annotate

             infile_basename_suffix_CHD_singles=.gridss.CHD_955_genes.vpot_final_output_file.vpot_gt_0.txt
outfile_basename_suffix_CHD_singles=_"${batch}"__gridss.CHD_955_genes.vpot_final_output_file.vpot_gt_0.txt

             infile_basename_suffix_CHD_family=.gridss.CHD_955_genes.vpot_final_output_file.txt
outfile_basename_suffix_CHD_family=_"${batch}"__gridss.CHD_genes.vpot_final_output_file.txt

             infile_basename_suffix_DCM_singles=.gridss.DCM_909_genes.vpot_final_output_file.vpot_gt_0.txt
outfile_basename_suffix_DCM_singles=_"${batch}"__gridss.DCM_909_genes.vpot_final_output_file.vpot_gt_0.txt

             infile_basename_suffix_DCM_family=.gridss.DCM_909_genes.vpot_final_output_file.txt
outfile_basename_suffix_DCM_family=_"${batch}"__gridss.DCM_909_genes.vpot_final_output_file.txt

tmpdir=/my/directory/variant_calling_and_annotation_pipeline/code/tmp_for_gridss_7_merge_samples
mkdir -p "${tmpdir}"

cohorts=(CHD DCM)

# Call vpot for singletons, and also for famiies.
# Merged singletons report will get vpot-prioritised variants with vpot score > 0.
# Merged families report will get all vpot-priroitised variants for manual inheritance modelling given that there are no many to go through.

for cohort in "${cohorts[@]}"; do

  echo 'Processing:' $cohort

  for singles_or_family in "singles" "family"; do

    echo '   Sub-processing:' $singles_or_family

    infile_basename_suffix_varname=infile_basename_suffix_"${cohort}"_"${singles_or_family}"
    outfile_basename_suffix_varname=outfile_basename_suffix_"${cohort}"_"${singles_or_family}"

    infile_basename_suffix="${!infile_basename_suffix_varname}"
    outfile_basename_suffix="${!outfile_basename_suffix_varname}"

    echo 'infile_basename_suffix_varname' $infile_basename_suffix_varname
    echo ''
    echo 'outfile_basename_suffix_varname' $outfile_basename_suffix_varname
    echo ''
    echo 'infile_basename_suffix' $infile_basename_suffix
    echo ''
    echo 'outfile_basename_suffix' $outfile_basename_suffix
    echo ''

    outfile="${indir}"/"${cohort}_${singles_or_family}_${outfile_basename_suffix}"
    string_of_files=''

    if [[ $infile_basename_suffix != "" ]]; then

      # read through the manifest list of families or singles, each one will have a structural variant file,
      # then see if it is in the cohort and singles_or_family that we are currently processing to merge into one file
      one_file=''
      while IFS= read -r inline; do

        IFS=$'\t' read -r -a array <<< "$inline"
        this_flagship="${array[1]}"
        this_sample="${array[3]}"
        this_cohort="${array[0]}"
        this_family="${array[2]}"
        if [[ $this_family == "." ]]; then
          # this manifest entry is for a singleton
          this_singles_or_family="singles"
        else
          # this manifest entry is for a family
          this_singles_or_family="family"
        fi

        if [[ "$this_cohort" == "$cohort" ]] && [[ "$this_singles_or_family" == "$singles_or_family" ]]; then

          # This family/sample is in this cohort, and is of the type (family of singleton) that we are currently processing to merge into one file

          if [[ $this_singles_or_family == "family" ]]; then
            # get the structural variants file for this family
            infile=$(ls -1 "${indir}"/"${this_family}"_*"${infile_basename_suffix}" | head -n 1)
            if [[ -f "$infile" ]]; then
              echo 'Found trio or family file for' $this_family ':' $infile
            else
              echo 'ERROR: Cannot find trio or family file for' $this_family
            fi

          else # [[ $this_singles_or_family == "singles" ]]; then
            # get the structural variants file for this singleton
            infile=$(ls -1 "${indir}"/"${this_sample}""${infile_basename_suffix}" | head -n 1)
            if [[ -f "$infile" ]]; then
              echo 'Found singleton file for' $this_sample ':' $infile
            else
              echo "ERROR: Cannot find file for" $this_sample
            fi
          fi

          echo ''
          echo 'Input file that will be used:' $infile
          echo ''

          # If this cohort can have trios, then expand number of sample fields to 6 if not already 6 (to allow trios plus three more family members)
          # So there will be 6 sample GT fields, one INFO_VAF field, 6 sample VAF fields, which is 13 fields.
          if [[ $singles_or_family == "family" ]]; then

            infile_basename=$(basename $infile)
            tmpfile="${tmpdir}"/"${infile_basename}"

            cols_string=$(head -n 1 $infile)
            IFS=$'\t' read -r -a array <<< "$cols_string"
            i=0
            col_sample_id=-1
            col_before_samples=-1
            col_vaf_in_info=-1
            col_after_samples=-1
            for this_col in "${array[@]}"; do
              i=$(( i + 1 ))
              if [[ $this_col == "Sample_Id" ]]; then # we need to repeat this column of sample_ids next to the GT and VAF columns 
                col_sample_id=$i
              fi
              if [[ $this_col == "INVBND2OVERLAP" ]]; then # in gridss file, this is the column before sample GTs
                col_before_samples=$i
              fi
              if [[ $this_col == "INV5" ]]; then # in manta file, this is the column before sample GTs
                col_before_samples=$i
              fi
              if [[ $this_col == "VAF_in_INFO" ]]; then # in gridss and manta files, this is the column after sample GTs and before sample VAFs
                col_vaf_in_info=$i
              fi
              if [[ $this_col == "gene_CDSexons" ]]; then # in gridss and manta files, this is the column after sample VAFs
                col_after_samples=$i
              fi
            done
            num_samples=$(( col_vaf_in_info - col_before_samples - 1 ))
            echo 'infile' $infile
            echo 'col_sample_id' $col_sample_id
            echo 'col_before_samples' $col_before_samples
            echo 'col_vaf_in_info' $col_vaf_in_info
            echo 'col_after_samples' $col_after_samples
            echo 'num_samples' $num_samples

            echo 'awk -v max_num_samples=10 -v col_sample_id='"$col_sample_id" '-v num_samples='"$num_samples" '-v col_before_samples='"$col_before_samples" '-v col_vaf_in_info='"$col_vaf_in_info" '-v col_after_samples='"$col_after_samples" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {printf $1; for (i=2; i<=col_before_samples; i++) printf FS $i; printf FS "family_sample_ids" FS "sample1_genotype" FS "sample2_genotype" FS "sample3_genotype" FS "sample4_genotype" FS "sample5_genotype" FS "sample6_genotype" FS "sample7_genotype" FS "sample8_genotype" FS "sample9_genotype" FS "sample10_genotype" FS "all_genotypes"; printf FS $col_vaf_in_info; printf FS "sample1_vaf" FS "sample2_vaf" FS "sample3_vaf" FS "sample4_vaf" FS "sample5_vaf" FS "sample6_vaf" FS "sample7_vaf" FS "sample8_vaf" FS "sample9_vaf" FS "sample10_vaf"; for (i=col_after_samples; i<=NF; i++) printf FS $i; printf "\n"} else {printf $1; for (i=2; i<=col_before_samples; i++) {printf FS $i}; printf FS $col_sample_id; num_missing=max_num_samples-num_samples; all_genotypes=""; for (i=(col_before_samples+1); i<=(col_before_samples+num_samples); i++) {printf FS $i; all_genotypes=all_genotypes""$i"_"}; if (num_missing>0) {for (i=1; i<=num_missing; i++) {printf FS "."}; all_genotypes=all_genotypes".""_"}; printf FS all_genotypes; printf FS $col_vaf_in_info; for (i=(col_vaf_in_info+1); i<=(col_vaf_in_info+num_samples); i++) {printf FS $i}; if (num_missing>0) {for (i=1; i<=num_missing; i++) {printf FS "."}}; for (i=col_after_samples; i<=NF; i++) {printf FS $i}; printf "\n"}}' $infile '>' $tmpfile

            awk -v max_num_samples=10 -v col_sample_id="$col_sample_id" -v num_samples="$num_samples" -v col_before_samples="$col_before_samples" -v col_vaf_in_info="$col_vaf_in_info" -v col_after_samples="$col_after_samples" 'BEGIN {FS="\t";OFS="\t"} {if (NR==1) {printf $1; for (i=2; i<=col_before_samples; i++) printf FS $i; printf FS "family_sample_ids" FS "sample1_genotype" FS "sample2_genotype" FS "sample3_genotype" FS "sample4_genotype" FS "sample5_genotype" FS "sample6_genotype" FS "sample7_genotype" FS "sample8_genotype" FS "sample9_genotype" FS "sample10_genotype" FS "all_genotypes"; printf FS $col_vaf_in_info; printf FS "sample1_vaf" FS "sample2_vaf" FS "sample3_vaf" FS "sample4_vaf" FS "sample5_vaf" FS "sample6_vaf" FS "sample7_vaf" FS "sample8_vaf" FS "sample9_vaf" FS "sample10_vaf"; for (i=col_after_samples; i<=NF; i++) printf FS $i; printf "\n"} else {printf $1; for (i=2; i<=col_before_samples; i++) {printf FS $i}; printf FS $col_sample_id; num_missing=max_num_samples-num_samples; all_genotypes=""; for (i=(col_before_samples+1); i<=(col_before_samples+num_samples); i++) {printf FS $i; all_genotypes=all_genotypes""$i"_"}; if (num_missing>0) {for (i=1; i<=num_missing; i++) {printf FS "."}; all_genotypes=all_genotypes".""_"}; printf FS all_genotypes; printf FS $col_vaf_in_info; for (i=(col_vaf_in_info+1); i<=(col_vaf_in_info+num_samples); i++) {printf FS $i}; if (num_missing>0) {for (i=1; i<=num_missing; i++) {printf FS "."}}; for (i=col_after_samples; i<=NF; i++) {printf FS $i}; printf "\n"}}' $infile > $tmpfile

            infile=$tmpfile
          fi
          # else, we are processing singletons, so no need to change the structural variants file format of singletons

          new_string_of_files="${string_of_files} ${infile}"
          string_of_files=$new_string_of_files
          one_file=$infile

        fi
      done < "$in_manifest"

      if [[ $one_file != '' ]]; then # some files were found in this cohort and singles_or_family, so output a merged file of those files

        tmphdr="${tmpdir}"/gridss_7_merge_samples.temp_gridss_hdr.txt
        tmphdr2="${tmpdir}"/gridss_7_merge_samples.temp_gridss_hdr2.txt

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

