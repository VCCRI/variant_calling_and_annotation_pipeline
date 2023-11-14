# variant_calling_and_annotation_pipeline
Scripts to be used on a high performance computing (HPC) PBS system to call, annotate and analyse variants for multiple samples of whole genome sequencing (WGS).

To install the scripts and some of the data and components of this pipeline, follow the instructions in list_tasks_to_install_pipeline.txt

Please note that this pipeline consists of mainly scripts to manage the running and reporting of the various pipeline softwares.
This pipeline does not contain the softwares that the scripts run nor the reference datasets that those softwares and the scripts use.
Those softwares and reference data need to be obtained and installed separately, and the system location pointed to by the appropriate variables in the script code/where_are_softwares_and_references.sh
Softwares called by the pipeline script and that need to be obtained and installed separately include: bwa, GATK, platypus, annovar, vep, VPOT, gridss, manta, cnvnator, conanvarvar, and spliceogen.
Reference data that need to be obtained and installed separately include: human reference genome, annovar annotation file, vep annotation file, gnomad population data, and spliceai data.

To run this pipeline to align Illumina paired-end sequencing reads to the human genome and then call and annotate variants, follow the instructions in list_jobs_to_run_pipeline.txt
