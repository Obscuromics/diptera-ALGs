#!/bin/bash

#BSUB -o logs/infer_alg.out.%J
#BSUB -e logs/infer_alg.err.%J
#BSUB -q week
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

module load nextflow/24.04.2-5914
module load mafft/7.525--h031d066_1
module load busco/5.8.2--pyhdfd78af_0
module load trimal/1.4.1--h4ac6f70_8
module load syngraph/0.0.1-c3

nextflow run $NF_PATH/infer_alg/infer_alg_main.nf -params-file $NF_PATH/infer_alg/diptera_alg_params.json -c $NF_PATH/conf/nextflow.config -with-conda /software/treeoflife/conda/users/envs/team360/se13/infer_alg -w data/workdir/infer_alg -resume  # -with-report -with-trace -with-timeline -with-dag
