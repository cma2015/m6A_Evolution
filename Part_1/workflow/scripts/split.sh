set -ue

# Get the bam file from the command line
DATA=$1
OUTPATH=$2
# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
## -f contain -F remove  128: second in pair 16: read reverse strand
## 80: read reverse strand  first in pair
## 144:  read reverse strand second in pair
## 64: first in pair
samtools view -@ 5 -b -f 128 -F 16 $DATA > ${OUTPATH}_fwd1.bam &&
samtools index -@ 5 ${OUTPATH}_fwd1.bam &&

samtools view -@ 5 -b -f 80 $DATA > ${OUTPATH}_fwd2.bam &&
samtools index -@ 5 ${OUTPATH}_fwd2.bam &&

#
# Combine alignments that originate on the forward strand.
#
samtools merge -@ 5 -f ${OUTPATH}_fwd.bam ${OUTPATH}_fwd1.bam ${OUTPATH}_fwd2.bam &&
samtools index -@ 5 ${OUTPATH}_fwd.bam &&

rm ${OUTPATH}_fwd1.bam ${OUTPATH}_fwd2.bam ${OUTPATH}_fwd1.bam.bai ${OUTPATH}_fwd2.bam.bai &
# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -@ 5 -b -f 144 $DATA > ${OUTPATH}_rev1.bam &&
samtools index -@ 5 ${OUTPATH}_rev1.bam &&

samtools view -@ 5 -b -f 64 -F 16 $DATA > ${OUTPATH}_rev2.bam &&
samtools index -@ 5 ${OUTPATH}_rev2.bam &&
#
# Combine alignments that originate on the reverse strand.
#
samtools merge -@ 5 -f ${OUTPATH}_rev.bam ${OUTPATH}_rev1.bam ${OUTPATH}_rev2.bam &&
samtools index -@ 5 ${OUTPATH}_rev.bam &&

rm ${OUTPATH}_rev1.bam ${OUTPATH}_rev2.bam ${OUTPATH}_rev1.bam.bai ${OUTPATH}_rev2.bam.bai &

wait