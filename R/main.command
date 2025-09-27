#!/bin/bash

bash scripts/OT_detector.sh
bash scripts/OT_mapper.sh analysis/OT/OT_candidate.bed
Rscript scripts/primer_generate.R -i analysis/OT_mapper_results/OT_mapped.tsv -o analysis/OT_mapper_results/OT_with_primer.tsv
