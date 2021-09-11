# Time: 20200706
# athCrBgd-m6A-seq analysis pipline with snakemake
SPECIES = ['athCrBgd1', 'athCrBgd2']
TREAT = ['IP', 'Input']
REP = [1, 2]

rule all:
    input:
        expand('{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam',species = SPECIES, treat = TREAT, rep = REP),
        expand('{species}/03BED/{species}_{treat}{rep}.bed',species = SPECIES, treat = TREAT, rep = REP),
        expand('{species}/04Peaks/{species}_peaks_R{rep}.txt',species=SPECIES,treat=TREAT,rep = REP),
        expand('{species}/04Peaks/{species}_peaks_R{rep}.bed',species = SPECIES,treat = TREAT,rep = REP),
        expand("{species}/04Peaks/{species}_peaksGene_tmp_R{rep}.bed", species = SPECIES, rep = REP),
        expand("{species}/04Gene/{species}_{treat}{rep}_gene.txt", species = SPECIES, treat = TREAT, rep = REP),
        expand("{species}_corrplot.pdf", species = SPECIES),
        expand("{species}/05Ratio/{species}_gene_ratio.txt", species = SPECIES)


rule star:
    input:
        reads = '{species}/01cleandata/{species}_{treat}{rep}.fastq.gz',
    output:
        bam = '{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam',
    threads: 2
    shell:
        """
        set -x
        STAR    --runThreadN 15 \
                --genomeDir ../genome/ath/STAR_index_100 \
                --readFilesIn {input} \
                --readFilesCommand zcat \
                --bamRemoveDuplicatesType UniqueIdentical \
                --limitBAMsortRAM 3762080671 \
                --outFileNamePrefix {wildcards.species}/02STAR/{wildcards.species}_{wildcards.treat}{wildcards.rep}_ \
                --outReadsUnmapped Fastx \
                --outSAMtype BAM   SortedByCoordinate  \
                --outSAMattrIHstart 0 \
                --outFilterMultimapNmax 5 \
                --outFilterMismatchNmax 0 \
                --clip3pAdapterSeq TGGAATTCTCGGGTGCCAAGG \
                --alignIntronMin 10 \
                --alignIntronMax 2000
        """

rule bamtobed:
    input:
        bamfile = '{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam',
        gfffile = '../genome/ath/ath.gtf'
    output:'{species}/03BED/{species}_{treat}{rep}.bed'
    threads: 1
    shell:
        """
        set -x
        file_prefix={wildcards.species}/03BED/{wildcards.species}_{wildcards.treat}{wildcards.rep}
        samtools view -F 1796 -q 30 -u {input.bamfile} | samtools sort -n - -o ${{file_prefix}}_final.bam
        bedtools bamtobed -i ${{file_prefix}}_final.bam > ${{file_prefix}}_tmp.bed
        bedtools sort -i {input.gfffile} | bedtools merge -i - >${{file_prefix}}_gene.bed
        bedtools intersect -wa -a ${{file_prefix}}_tmp.bed -b ${{file_prefix}}_gene.bed -f 0.5 | sort -k 1,1 -k 2,2n | uniq > {output}
        """

rule peakcalling:
    input:
        IP = "{species}/03BED/{species}_IP{rep}.bed",
        Input = "{species}/03BED/{species}_Input{rep}.bed",
        genome = "../genome/ath/genome.size"
    output: "{species}/04Peaks/{species}_peaks_R{rep}.txt"
    threads: 1
    shell:
        """
        set -x
        which Rscript
        /usr/bin/Rscript ../scripts/PEA_peak_calling.R -i {input.IP} -n {input.Input} -g {input.genome} -o {output}
        """

rule peak2bed:
    input: "{species}/04Peaks/{species}_peaks_R{rep}.txt"
    output: "{species}/04Peaks/{species}_peaks_R{rep}.bed"
    threads: 1
    shell:
        """
        awk '{{printf "%s\\t%d\\t%d\\n",$1,$2-1,$3}}' {input} > {output}
        """

rule stringtie:
    input:
        bam = "{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam",
        gtf = "../genome/ath/ath.gtf"
    output:
        gtf_s = "{species}/04Gene/{species}_{treat}{rep}.gtf",
        gene_exp = "{species}/04Gene/{species}_{treat}{rep}_gene.txt"
    threads: 2
    shell: "stringtie -e -p 8 -G {input.gtf} -j 10 -o {output.gtf_s} -A {output.gene_exp} {input.bam}"

rule exon_gff3:
    input: "../genome/ath/ath.gff3"
    output: "../genome/ath/ath.exons.gff3"
    shell: "awk '$3==\"exon\"||$3==\"CDS\"' {input} | sort -k1,1 -k4,4n >{output}"


rule add_gene_strand:
    input:
        peak = "{species}/04Peaks/{species}_peaks_R{rep}.bed",
        gene = "../genome/ath/ath.bed",
        expGene = "{species}/04Gene/{species}_Input{rep}_gene.txt",
        exongff = "../genome/ath/ath.exons.gff3"
    output:
        geneexp = "{species}/04Gene/{species}_Input{rep}_feature.txt",
        genepeak = "{species}/04Peaks/{species}_peaksGene_R{rep}.txt"
    threads: 1
    shell:
        """
        awk -F"\t" 'NR==FNR{{a[$1]=$9;next}}{{OFS="\t";if(a[$4]){{print $0"\t"a[$4]}}else{{print $0"\t"0}}}}' {input.expGene} {input.gene} >{output.geneexp}
        bedtools intersect -wao -a {input.peak} -b {output.geneexp} | cut -f 1,2,3,7-11 >{output.genepeak}
        """

rule remove_overlap_gene:
    input: "{species}/04Peaks/{species}_peaksGene_R{rep}.txt"
    output:
        bed = "{species}/04Peaks/{species}_peaksGene_R{rep}.bed"
    threads: 1
    shell: "/usr/bin/Rscript ../scripts/peak_remove.R {input} 1 {output.bed}"

rule remove_cut:
    input: "{species}/04Peaks/{species}_peaksGene_R{rep}.bed"
    output: "{species}/04Peaks/{species}_peaksGene_tmp_R{rep}.bed"
    shell: "cut -f 1-3 {input} | sort -k1,1 -k2,2n | uniq > {output}"

rule diffbind_input:
    input: "{species}/04Peaks/{species}_peaksGene_R1.bed",
            "{species}/04Peaks/{species}_peaksGene_R2.bed"
    output: "{species}/DiffBindPath.csv"
    shell: "sed 's/argPath/{wildcards.species}/g;s/data_12ss/data_public/g' ../scripts/DiffBindPathP2.csv > {wildcards.species}/DiffBindPath.csv"

rule corrplot:
    input: "{species}/DiffBindPath.csv"
    output:
        vennplot = "{species}_venn.pdf",
        corrplot = "{species}_corrplot.pdf"
    shell: "/usr/bin/Rscript ../scripts/diffbind.R {input} {wildcards.species} {output.vennplot} {output.corrplot}"

rule intersect_peak:
    input: "{species}/04Peaks/{species}_peaksGene_tmp_R1.bed",
            "{species}/04Peaks/{species}_peaksGene_tmp_R2.bed"
    output: "{species}/04Peaks/{species}_peaksGene_tmp.bed"
    threads: 1
    shell:
        """
        /usr/bin/Rscript ../scripts/samples_merge.R {input} {output}
        """

rule gene_expression:
    input: "{species}/04Gene/{species}_Input1_gene.txt",
            "{species}/04Gene/{species}_Input2_gene.txt"
    output: "{species}/04Gene/{species}_gene_expression.txt"
    threads: 1
    shell:
        """
        /usr/bin/Rscript ../scripts/gene_replication_merge.R {wildcards.species}/04Gene {wildcards.species}_Input {output}
        """

rule merge_gene:
    input:
        peak = "{species}/04Peaks/{species}_peaksGene_tmp.bed",
        gene = "../genome/ath/ath.bed",
        expGene = "{species}/04Gene/{species}_gene_expression.txt"
    output:
        geneexp = "{species}/04Gene/{species}_gene_feature.txt",
        overlap = "{species}/04Peaks/{species}_merge.txt"
    threads: 1
    shell:
        """
        awk -F"\t" 'NR==FNR{{a[$1]=$NF;next}}{{OFS="\t";if(a[$4]){{print $0"\t"a[$4]}}else{{print $0"\t"0}}}}' {input.expGene} {input.gene} >{output.geneexp}
        sort -k1,1 -k2,2n {output.geneexp} | bedtools intersect -wao -a {input.peak} -b - | cut -f 1,2,3,10-14 >{output.overlap}
        """

rule merge_remove:
    input: "{species}/04Peaks/{species}_merge.txt"
    output:
        bed = "{species}/04Peaks/{species}_merge.bed"
    threads: 1
    shell: "/usr/bin/Rscript ../scripts/peak_remove2.R {input} 1 {output.bed}"


rule peak_cut:
    input: "{species}/04Peaks/{species}_merge.bed"
    output: "{species}/04Peaks/{species}_tmp_merge.bed"
    shell: "cut -f 1-3 {input} > {output}"


rule peak_RPKM:
    input:
        dbPath = "{species}/DiffBindPath.csv",
        mergePeak = "{species}/04Peaks/{species}_tmp_merge.bed"
    output: "{species}/05Ratio/{species}_peak_ratio.txt"
    shell: "/usr/bin/Rscript ../scripts/diffbind2.R {input.dbPath} {input.mergePeak} {output}"


rule gene_RPKM:
    input:
        peakRatio = "{species}/05Ratio/{species}_peak_ratio.txt",
        mergePeak = "{species}/04Peaks/{species}_merge.bed"
    output: "{species}/05Ratio/{species}_gene_ratio.txt"
    shell:
        """
        awk -F"\t" 'NR==FNR{{a[$1"_"$2"_"$3]=$0;next}}{{print($0"\t"a[$1"_"$2"_"$3])}}' {input.peakRatio} {input.mergePeak} | \
        cut -f 1-6,10-15 | awk '{{OFS="\t";$2+=1;print $0}}' - >{output}
        """
