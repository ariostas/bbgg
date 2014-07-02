#!/bin/bash

echo "Submiting all jobs"

bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"Signal\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"B\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"BB\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"BBB\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"Bj\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"Bjj\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"H\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"LL\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"LLB\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"tB\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"tj\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"tt\"\)
bsub -q 1nd -W 1250 run.sh bbgg.C+\(\"ttB\"\)

