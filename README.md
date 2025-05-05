# variant_calling_and_annotation_pipeline
Scripts to be used on a high performance computing (HPC) PBS system to call, annotate and analyse variants for multiple samples of whole genome sequencing (WGS).

To install the scripts and some of the data and components of this pipeline, follow the instructions in list_tasks_to_install_pipeline.txt

Please note that this pipeline consists of mainly scripts to manage the running and reporting of the various pipeline softwares.
This pipeline does not contain the softwares that the scripts run nor the reference datasets that those softwares and the scripts use.
Those softwares and reference data need to be obtained and installed separately, and the system location pointed to by the appropriate variables in the script code/where_are_softwares_and_references.sh
Softwares called by the pipeline script and that need to be obtained and installed separately include: bwa, GATK, platypus, annovar, vep, VPOT, gridss, manta, cnvnator, conanvarvar, and spliceogen.
Reference data that need to be obtained and installed separately include: human reference genome, annovar annotation file, vep annotation file, gnomad population data, and spliceai data.
For the list of softwares that need to be installed, please look at code/where_are_softwares_and_references.sh. Every variable that points to a software path is a software that needs to be installed.

To run this pipeline to align Illumina paired-end sequencing reads to the human genome and then call and annotate variants, follow the instructions in list_jobs_to_run_pipeline.txt

Notes about VPOT and pedigree file:
This pipeline runs the VPOT software () to output variants that conforms to various types of inheritance models.<br>
Eg. If the child is affected by the genetic disease and the parents are not, then run VPOT for the following inheritance models: DN (de-novo), CH (compound heterozygous), AR (autosomal recessive).<br>
Eg. If the child and mother are affected and the father is not, then run VPOT for: AD-mother (autosomal dominant on the mother's side).<br>
Eg. If you have several family members and some are affected and others are not, then run Case-Control.<br>
VPOT uses the standard pedigree file to know who is affected and how they are related.<br>
Eg. Pedigree file for 2 families, each having affect child (one female, one male) and unaffected parents:<br>
<pre>
$ cat MY_SAMPLES.ped
#Family	Subject	Father	Mother	Sex	Phenotype
FAM1	SAMP1	SAMP2	SAMP3	2	2
FAM1	SAMP3	0	0	2	1
FAM1	SAMP2	0	0	1	1
FAM2	SAMP21	SAMP22	SAMP23	1	2
FAM2	SAMP22	0	0	1	1
FAM2	SAMP23	0	0	2	1
</pre>
For a duo consisting of an affected child and unaffected mother, one can run VPOT for DN and CH. Here is how to create the VPOT pedigree file:
<pre>
$ cat MY_SAMPLES3.ped
#Family Subject Father  Mother  Sex     Phenotype
FAM3	SAMP31	0	SAMP32	2	2
FAM3	SAMP32	0	0	2	1
FAM3	0	0	0	1	1
</pre>

