#!/bin/bash

# set environment variables for softwares and references
currdir=$(pwd)
sw_and_refs="${currdir}"/where_are_softwares_and_references.sh
. "${sw_and_refs}"
# export spliceogen_sw="${sw}"/spliceogen_2020july/

cp -r "${spliceogen_sw}"/Spliceogen/* ../working_directory/spliceogen
cp "${spliceogen_sw}"/resources/spliceogen_run_when_launched_in_pbs.sh ../working_directory/spliceogen

cp "${spliceogen_sw}"/resources/spliceogen_launch_one_sample_to_pbs.pbs .
cp "${spliceogen_sw}"/resources/spliceogen_run_to_launch_multiple_samples_to_pbs.sh .
cp "${spliceogen_sw}"/resources/spliceogen_run_when_launched_in_pbs.sh .

