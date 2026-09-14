#!/bin/bash

snakemake \
   -s Snakefile_mvmr_cf \
   --keep-going \
   --notemp \
   --jobs 96 \
   --max-jobs-per-second 5 \
   --latency-wait 30 \
   --default-resources mem_mb=5000 runtime=120 account="jvmorr0"\
   --executor slurm
