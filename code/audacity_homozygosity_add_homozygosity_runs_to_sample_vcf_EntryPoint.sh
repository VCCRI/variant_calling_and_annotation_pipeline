#!/bin/bash

manifest=$1

# NCI queue prioritises running jobs in the queue when the number of jobs queued is limited.
# So queue only a limited number of jobs, and keep interrogating the number of jobs queued so that when queued jobs start running (and thus are no longer queued), more jobs will be queued.
while true; do
  ./audacity_homozygosity_add_homozygosity_runs_to_sample_vcf_SubmitJobs.sh $manifest
  sleep 900
done


