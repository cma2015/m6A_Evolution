# Time: 2019-10
# m6A-seq analysis pipeline with snakemake

## Organize files for processing
FOLDER = ["."]
SPECIES = "maizePpHy"
TREAT = ["IP", "Input"]
REP = [1, 2]

rule fastq2peak:
    input:
        expand("{datadir}/{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam", datadir=FOLDER, species = SPECIES, treat = TREAT, rep = REP),
        expand("{datadir}/{species}/01multiQC", datadir=FOLDER, species = SPECIES),
        expand("{datadir}/{species}/04Peaks/{species}_peaksGene_tmp_R{rep}.bed", datadir=FOLDER, species = SPECIES, rep = REP),
        expand("{datadir}/{species}/04Gene/{species}_{treat}{rep}_gene.txt", datadir=FOLDER, species = SPECIES, treat = TREAT, rep = REP),
        expand("{datadir}/{species}_corrplot.pdf", datadir=FOLDER, species = SPECIES),
        expand("{datadir}/{species}/05Ratio/{species}_gene_ratio.txt", datadir=FOLDER, species = SPECIES)

rule fastp:
    input:
        reads1 = "{datadir}/{species}/00rawdata/{species}_{treat}{rep}_R1.fastq.gz",
        reads2 = "{datadir}/{species}/00rawdata/{species}_{treat}{rep}_R2.fastq.gz"
    output:
        creads1 = "{datadir}/{species}/01cleandata/{species}_{treat}{rep}_c_R1.fastq.gz",
        creads2 = "{datadir}/{species}/01cleandata/{species}_{treat}{rep}_c_R2.fastq.gz",
        json = "{datadir}/{species}/01fastp_html/{species}_{treat}{rep}_fastp.json",
        html = "{datadir}/{species}/01fastp_html/{species}_{treat}{rep}_fastp.html"
    threads: 1
    shell:
        """
        mkdir -p {wildcards.datadir}/{wildcards.species}/01fastqc
        fastp -w 10 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
        -i {input.reads1} -o {output.creads1} -I {input.reads2} -O {output.creads2} -j {output.json} -h {output.html}
        ~/publicdir/software/FastQC/fastqc -t 10 -O {wildcards.datadir}/{wildcards.species}/01fastqc {input.reads1}
        ~/publicdir/software/FastQC/fastqc -t 10 -O {wildcards.datadir}/{wildcards.species}/01fastqc {input.reads2}
        ~/publicdir/software/FastQC/fastqc -t 10 -O {wildcards.datadir}/{wildcards.species}/01fastqc {output.creads1}
        ~/publicdir/software/FastQC/fastqc -t 10 -O {wildcards.datadir}/{wildcards.species}/01fastqc {output.creads2}
        """

rule multiQC:
    output: "{datadir}/{species}/01multiQC"
    log: "{datadir}/{species}/01multiQC/multiqc.log"
    message: "multiqc for all logs"
    threads: 1
    shell:
        """
        multiqc {wildcards.datadir}/{wildcards.species}/01fastqc -o {output} -d -f -v -n multiQC 2> {log}
        """

rule star:
    input:
        reads1 = "{datadir}/{species}/01cleandata/{species}_{treat}{rep}_c_R1.fastq.gz",
        reads2 = "{datadir}/{species}/01cleandata/{species}_{treat}{rep}_c_R2.fastq.gz"
    output:
        bam = "{datadir}/{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam",
        align = "{datadir}/{species}/02STAR/{species}_{treat}{rep}_Log.final.out"
    params:
        index="../genome/zma/STAR_index"
    threads: 2
    shell:
        """
        set -x
        STAR --genomeDir {params.index} --runThreadN 15 --readFilesIn {input.reads1} {input.reads2} \
            --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
            --alignIntronMin 10 --alignIntronMax 17000 \
            --outFileNamePrefix {wildcards.datadir}/{wildcards.species}/02STAR/{wildcards.species}_{wildcards.treat}{wildcards.rep}_ \
            --outFilterMismatchNmax 0 --outSAMattrIHstart 0 --outFilterMultimapNmax 1 --outReadsUnmapped Fastx
        samtools index {output.bam}
        """

rule genome_fai:
    input: "../genome/zma/zma.fa"
    output: "../genome/zma/genome.size"
    threads: 1
    shell:
        """
        samtools faidx {input}
        cut -f 1,2 {input}.fai > {output}
        """

rule bamtobed:
    input:
        bamfile = "{datadir}/{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam",
        gfffile = "../genome/zma/zma.gtf"
    output: "{datadir}/{species}/03BED/{species}_{treat}{rep}.bed"
    threads: 1
    shell:
        """
        file_prefix={wildcards.datadir}/{wildcards.species}/03BED/{wildcards.species}_{wildcards.treat}{wildcards.rep}
        samtools view -F 1804 -f 2 -q 30 -u {input.bamfile} | samtools sort -n - -o ${{file_prefix}}_sort.bam
        samtools fixmate -r ${{file_prefix}}_sort.bam ${{file_prefix}}.fixmate.tmp
        samtools view -F 1804 -f 2 -u ${{file_prefix}}.fixmate.tmp | samtools sort -n - -o ${{file_prefix}}_final.bam
        bedtools bamtobed -bedpe -mate1 -i ${{file_prefix}}_final.bam | \
        awk '{{a=i++;if($2>=$5){{min=$5}}else{{min=$2}};if($3>$6){{max=$3}}else{{max=$6}};OFS="\t";print $1,min,max,"read_"a+1}}' - > ${{file_prefix}}_tmp.bed
        bedtools sort -i {input.gfffile} | bedtools merge -i - >${{file_prefix}}_gene.bed
        bedtools intersect -wa -a ${{file_prefix}}_tmp.bed -b ${{file_prefix}}_gene.bed -f 0.5 | sort -k 1,1 -k 2,2n | uniq > {output}
        """

rule peakcalling:
    input:
        IP = "{datadir}/{species}/03BED/{species}_IP{rep}.bed",
        Input = "{datadir}/{species}/03BED/{species}_Input{rep}.bed",
        genome = "../genome/zma/genome.size"
    output: "{datadir}/{species}/04Peaks/{species}_peaks_R{rep}.txt"
    threads: 1
    shell:
        """
        /usr/bin/Rscript ../scripts/PEA_peak_calling.R -i {input.IP} -n {input.Input} -g {input.genome} -o {output}
        """

rule peak2bed:
    input: "{datadir}/{species}/04Peaks/{species}_peaks_R{rep}.txt"
    output: "{datadir}/{species}/04Peaks/{species}_peaks_R{rep}.bed"
    threads: 1
    shell:
        """
        awk '{{printf "%s\\t%d\\t%d\\n",$1,$2-1,$3}}' {input} > {output}
        """

rule stringtie:
    input:
        bam = "{datadir}/{species}/02STAR/{species}_{treat}{rep}_Aligned.sortedByCoord.out.bam",
        gtf = "../genome/zma/zma.gtf"
    output:
        gtf_s = "{datadir}/{species}/04Gene/{species}_{treat}{rep}.gtf",
        gene_exp = "{datadir}/{species}/04Gene/{species}_{treat}{rep}_gene.txt"
    threads: 2
    shell: "stringtie -e -p 8 -G {input.gtf} -j 10 -o {output.gtf_s} -A {output.gene_exp} {input.bam}"

rule exon_gff3:
    input: "../genome/zma/zma.gff3"
    output: "../genome/zma/zma.exons.gff3"
    shell: "awk '$3==\"exon\"||$3==\"CDS\"' {input} | sort -k1,1 -k4,4n >{output}"

rule add_gene_strand:
    input:
        peak = "{datadir}/{species}/04Peaks/{species}_peaks_R{rep}.bed",
        gene = "../genome/zma/zma.bed",
        expGene = "{datadir}/{species}/04Gene/{species}_Input{rep}_gene.txt",
        exongff = "../genome/zma/zma.exons.gff3"
    output:
        geneexp = "{datadir}/{species}/04Gene/{species}_Input{rep}_feature.txt",
        genepeak = "{datadir}/{species}/04Peaks/{species}_peaksGene_R{rep}.txt"
    threads: 1
    shell:
        """
        awk -F"\t" 'NR==FNR{{a[$1]=$9;next}}{{OFS="\t";if(a[$4]){{print $0"\t"a[$4]}}else{{print $0"\t"0}}}}' {input.expGene} {input.gene} >{output.geneexp}
        bedtools intersect -wao -a {input.peak} -b {output.geneexp} | cut -f 1,2,3,7-11 >{output.genepeak}
        """

rule remove_overlap_gene:
    input: "{datadir}/{species}/04Peaks/{species}_peaksGene_R{rep}.txt"
    output:
        bed = "{datadir}/{species}/04Peaks/{species}_peaksGene_R{rep}.bed"
    threads: 1
    shell: "/usr/bin/Rscript ../scripts/peak_remove.R {input} 1 {output.bed}"

rule remove_cut:
    input: "{datadir}/{species}/04Peaks/{species}_peaksGene_R{rep}.bed"
    output: "{datadir}/{species}/04Peaks/{species}_peaksGene_tmp_R{rep}.bed"
    shell: "cut -f 1-3 {input} | sort -k1,1 -k2,2n | uniq > {output}"

rule diffbind_input:
    input: "{datadir}/{species}/04Peaks/{species}_peaksGene_R1.bed",
            "{datadir}/{species}/04Peaks/{species}_peaksGene_R2.bed"
    output: "{datadir}/{species}/DiffBindPath.csv"
    shell: "sed 's/argPath/{wildcards.species}/g' ../scripts/DiffBindPathP2.csv > {wildcards.datadir}/{wildcards.species}/DiffBindPath.csv"

rule corrplot:
    input: "{datadir}/{species}/DiffBindPath.csv"
    output:
        vennplot = "{datadir}/{species}_venn.pdf",
        corrplot = "{datadir}/{species}_corrplot.pdf"
    shell: "/usr/bin/Rscript ../scripts/diffbind.R {input} {wildcards.species} {output.vennplot} {output.corrplot}"

rule intersect_peak:
    input: "{datadir}/{species}/04Peaks/{species}_peaksGene_tmp_R1.bed",
            "{datadir}/{species}/04Peaks/{species}_peaksGene_tmp_R2.bed"
    output: "{datadir}/{species}/04Peaks/{species}_peaksGene_tmp.bed"
    threads: 1
    shell:
        """
        /usr/bin/Rscript ../scripts/samples_merge.R {input} {output}
        """

rule gene_expression:
    input: "{datadir}/{species}/04Gene/{species}_Input1_gene.txt",
            "{datadir}/{species}/04Gene/{species}_Input2_gene.txt"
    output: "{datadir}/{species}/04Gene/{species}_gene_expression.txt"
    threads: 1
    shell:
        """
        /usr/bin/Rscript ../scripts/gene_replication_merge.R {wildcards.datadir}/{wildcards.species}/04Gene {wildcards.species}_Input {output}
        """

rule merge_gene:
    input:
        peak = "{datadir}/{species}/04Peaks/{species}_peaksGene_tmp.bed",
        gene = "../genome/zma/zma.bed",
        expGene = "{datadir}/{species}/04Gene/{species}_gene_expression.txt"
    output:
        geneexp = "{datadir}/{species}/04Gene/{species}_gene_feature.txt",
        overlap = "{datadir}/{species}/04Peaks/{species}_merge.txt"
    threads: 1
    shell:
        """
        awk -F"\t" 'NR==FNR{{a[$1]=$NF;next}}{{OFS="\t";if(a[$4]){{print $0"\t"a[$4]}}else{{print $0"\t"0}}}}' {input.expGene} {input.gene} >{output.geneexp}
        sort -k1,1 -k2,2n {output.geneexp} | bedtools intersect -wao -a {input.peak} -b - | cut -f 1,2,3,10-14 >{output.overlap}
        """

rule merge_remove:
    input: "{datadir}/{species}/04Peaks/{species}_merge.txt"
    output:
        bed = "{datadir}/{species}/04Peaks/{species}_merge.bed"
    threads: 1
    shell: "/usr/bin/Rscript ../scripts/peak_remove2.R {input} 1 {output.bed}"


rule peak_cut:
    input: "{datadir}/{species}/04Peaks/{species}_merge.bed"
    output: "{datadir}/{species}/04Peaks/{species}_tmp_merge.bed"
    shell: "cut -f 1-3 {input} > {output}"


rule peak_RPKM:
    input:
        dbPath = "{datadir}/{species}/DiffBindPath.csv",
        mergePeak = "{datadir}/{species}/04Peaks/{species}_tmp_merge.bed"
    output: "{datadir}/{species}/05Ratio/{species}_peak_ratio.txt"
    shell: "/usr/bin/Rscript ../scripts/diffbind2.R {input.dbPath} {input.mergePeak} {output}"


rule gene_RPKM:
    input:
        peakRatio = "{datadir}/{species}/05Ratio/{species}_peak_ratio.txt",
        mergePeak = "{datadir}/{species}/04Peaks/{species}_merge.bed"
    output: "{datadir}/{species}/05Ratio/{species}_gene_ratio.txt"
    shell:
        """
        awk -F"\t" 'NR==FNR{{a[$1"_"$2"_"$3]=$0;next}}{{print($0"\t"a[$1"_"$2"_"$3])}}' {input.peakRatio} {input.mergePeak} | \
        cut -f 1-6,10-15 | awk '{{OFS="\t";$2+=1;print $0}}' - >{output}
        """
