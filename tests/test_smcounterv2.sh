#!/bin/bash
set -e
codedir=$1

python $codedir/run.py --runPath /home/qiauser/test_v2/ --bamFile /home/qiauser/test_v2/NB956-240-3-10_S1.highconfidence.bam \
       --bedTarget /home/qiauser/test_v2/high.confidence.variants.bed --outPrefix NB956-240-3-10_S1.test --nCPU 2 --minBQ 25 \
       --minMQ 50 --hpLen 8 --mismatchThr 6.0 --primerDist 2 --consThr 0.8 --rpu 7.6 --primerSide 1 --minAltUMI 3 --maxAltAllele 2 \
       --refGenome /srv/qgen/data/genome/hg19/ucsc.hg19.fa --repBed /srv/qgen/data/annotation/simpleRepeat.full.bed \
       --srBed /srv/qgen/data/annotation/SR_LC_SL.full.bed

python $codedir/tests/compare_outlong.py /home/qiauser/test_v2/NB956-240-3-10_S1.highconfidence.VariantList.long.txt \
       /home/qiauser/test_v2/intermediate/NB956-240-3-10_S1.test.VariantList.long.txt True True
