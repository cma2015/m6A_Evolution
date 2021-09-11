#!/bin/bash
processing of strand-specific MeRIP-seq/m6A-seq data.

set -x

for SP in ath gar ghi pvu gma sbi zma ata tdi tae osa ppa maizePgLyj1 maizePgLyj2 maizePgLyj3
do
	if [[ ${SP} == "maizePgLyj1" ]] || [[ ${SP} == "maizePgLyj2" ]] || [[ ${SP} == "maizePgLyj3" ]];then
		RAWSP=Zma
		DATAPATH=data_public
		SAPNUM=2
	else
		RAWSP=${SP^}
		DATAPATH=data_12ss
		SAPNUM=3
	fi

    mkdir -p data/${SP}
    mkdir -p data/${SP}/02BAMsplit/ data/${SP}/03Peak/
    mkdir -p data/${SP}/03Peak/raw data/${SP}/03Peak/bed

    BEDPATH=data/${SP}/03Peak/bed

    for rep in `seq ${SAPNUM}`
    do
		if [[ ${SP} == "ata" ]] || [[ ${SP} == "tdi" ]] || [[ ${SP} == "tae" ]];then
        sh scripts/split_c.sh ${DATAPATH}/${SP}/02STAR/${SP}_Input${rep}_Aligned.sortedByCoord.out.bam  data/${SP}/02BAMsplit/${SP}_Input${rep} &
        sh scripts/split_c.sh ${DATAPATH}/${SP}/02STAR/${SP}_IP${rep}_Aligned.sortedByCoord.out.bam  data/${SP}/02BAMsplit/${SP}_IP${rep} &
        wait
		else
        sh scripts/split.sh ${DATAPATH}/${SP}/02STAR/${SP}_Input${rep}_Aligned.sortedByCoord.out.bam  data/${SP}/02BAMsplit/${SP}_Input${rep} &
        sh scripts/split.sh ${DATAPATH}/${SP}/02STAR/${SP}_IP${rep}_Aligned.sortedByCoord.out.bam  data/${SP}/02BAMsplit/${SP}_IP${rep} &
        wait
		fi

      for ii in Input IP
      do
		for jj in fwd rev
		do
			DATAPREFIX=data/${SP}/02BAMsplit/${SP}_${ii}${rep}_${jj} &&
			samtools view -@ 10 -f 2 -F 1804 -q 30 -u ${DATAPREFIX}.bam | samtools sort -@ 10 -n - -o ${DATAPREFIX}_sort.bam &&
			samtools fixmate -r ${DATAPREFIX}_sort.bam ${DATAPREFIX}.fixmate.tmp &&
			samtools view -@ 10 -f 2 -F 1804 -u ${DATAPREFIX}.fixmate.tmp | samtools sort -@ 10 -n - -o ${DATAPREFIX}_final.bam &&
			bedtools bamtobed -bedpe -mate1 -i ${DATAPREFIX}_final.bam | \
			awk '{{a=i++;if($2>=$5){{min=$5}}else{{min=$2}};if($3>$6){{max=$3}}else{{max=$6}};OFS="\t";print $1,min,max,"read_"a+1}}' - | sort -k 1,1 -k 2,2n - > ${DATAPREFIX}_tmp.bed &
		done
      done

      wait

        ip_fwd=`wc -l data/${SP}/02BAMsplit/${SP}_IP${rep}_fwd_tmp.bed | cut -f 1 -d " "`
        ip_rev=`wc -l data/${SP}/02BAMsplit/${SP}_IP${rep}_rev_tmp.bed | cut -f 1 -d " "`
        input_fwd=`wc -l data/${SP}/02BAMsplit/${SP}_Input${rep}_fwd_tmp.bed | cut -f 1 -d " "`
        input_rev=`wc -l data/${SP}/02BAMsplit/${SP}_Input${rep}_rev_tmp.bed | cut -f 1 -d " "`
        ip_count=`expr ${ip_fwd} + ${ip_rev}`
        input_count=`expr ${input_fwd} + ${input_rev}`

        for kk in fwd rev
        do
            /usr/bin/Rscript scripts/PEA_peak_calling.R \
            -i data/${SP}/02BAMsplit/${SP}_IP${rep}_${kk}_tmp.bed \
            -n data/${SP}/02BAMsplit/${SP}_Input${rep}_${kk}_tmp.bed \
            -g ~/a2z/metaPlants/workflow/genome/${RAWSP}/Genome/${RAWSP}.genome.size \
            -p ${ip_count} -q ${input_count} \
            -o data/${SP}/03Peak/raw/${SP}_${kk}_R${rep}.txt &&
            awk -F"\t" '{OFS="\t";print $1,$2-1,$3}' data/${SP}/03Peak/raw/${SP}_${kk}_R${rep}.txt \
            > data/${SP}/03Peak/bed/${SP}_${kk}_R${rep}.bed &
        done

        wait
    done
done


for SP in  ath gar ghi pvu gma sbi zma ata tdi tae osa ppa # maizePgLyj1 maizePgLyj2 maizePgLyj3 #
do
	if [[ ${SP} == "maizePgLyj1" ]] || [[ ${SP} == "maizePgLyj2" ]] || [[ ${SP} == "maizePgLyj3" ]];then
		RAWSP=Zma
		DATAPATH=~/a2z/m6A_13spp/workflow/data_public
		SAPNUM=2
	else
		RAWSP=${SP^}
		DATAPATH=~/a2z/m6A_13spp/workflow/data_12ss
		SAPNUM=3
	fi

	PCGGENE=~/a2z/metaPlants/workflow/genome/${RAWSP}/Annotation/${RAWSP}
	mkdir -p ~/a2z/m6A_13spp/Part_1/workflow/data/${SP}/04Gene
	mkdir -p ~/a2z/m6A_13spp/Part_1/workflow/data/${SP}/05Filter
	BAMPATH=~/a2z/m6A_13spp/Part_1/workflow/data/${SP}/02BAMsplit
	PEAKPATH=~/a2z/m6A_13spp/Part_1/workflow/data/${SP}/03Peak
	GENEPATH=~/a2z/m6A_13spp/Part_1/workflow/data/${SP}/04Gene
	FILTERDATA=~/a2z/m6A_13spp/Part_1/workflow/data/${SP}/05Filter

	for rep in `seq ${SAPNUM}`
	do
		awk -F"\t" '{OFS="\t";print $0,".", ".", "+"}' \
		${PEAKPATH}/bed/${SP}_fwd_R${rep}.bed >${PEAKPATH}/${SP}_fwd_R${rep}.bed 
		awk -F"\t" '{OFS="\t";print $0,".", ".", "-"}' \
		${PEAKPATH}/bed/${SP}_rev_R${rep}.bed >${PEAKPATH}/${SP}_rev_R${rep}.bed 
	done

	for rep in `seq ${SAPNUM}`
	do
		DATAPREFIX=${DATAPATH}/${SP}/03BED/${SP}_Input${rep}_final
		samtools sort -@ 10 -o ${DATAPREFIX}_sort.bam ${DATAPREFIX}.bam
		stringtie -e -p 18 --rf \
		-G ${PCGGENE}.exons.gtf \
		-j 10 -o ${GENEPATH}/raw/${SP}_Input${rep}.gtf \
		-A ${GENEPATH}/raw/${SP}_Input${rep}_gene.txt ${DATAPREFIX}_sort.bam
		for strtype in fwd rev
		do
			awk -F"\t" 'NR>1&&$9>=1{OFS="\t";print $3,$5-1,$6,$1,$9,$4}' \
			${GENEPATH}/raw/${SP}_Input${rep}_gene.txt | \
			bedtools intersect -nonamecheck -wo -s -a ${PEAKPATH}/${SP}_${strtype}_R${rep}.bed -b - | \
			awk -F"\t" '$13>=100{OFS="\t";print $1,$2,$3,$10,$11,$12}' - | \
			awk -F"\t" 'NR==FNR{a[$4]=$4;next}{if(a[$4]){print $0}}' \
			${PCGGENE}.PCG.info.bed - > data/${SP}/05Filter/${SP}_${strtype}_R${rep}.bed
		done
	done

	for strtype in fwd rev
	do
	if [[ ${SAPNUM} == 3 ]];then
	/usr/bin/Rscript scripts/samples_peak_merge.R ${FILTERDATA}/${SP}_${strtype}_R1.bed ${FILTERDATA}/${SP}_${strtype}_R2.bed \
	${FILTERDATA}/${SP}_${strtype}_R3.bed ${FILTERDATA}/${SP}_${strtype}_merge.bed
	elif [[ ${SAPNUM} == 2 ]];then
	/usr/bin/Rscript scripts/samples_peak_merge.R ${FILTERDATA}/${SP}_${strtype}_R1.bed ${FILTERDATA}/${SP}_${strtype}_R2.bed ${FILTERDATA}/${SP}_${strtype}_merge.bed
	fi
	done

	awk -F"\t" '{OFS="\t";$6="+";print $0}' ${FILTERDATA}/${SP}_fwd_merge.bed >${FILTERDATA}/${SP}_fwd_merge_mod.bed
	awk -F"\t" '{OFS="\t";$6="-";print $0}' ${FILTERDATA}/${SP}_rev_merge.bed >${FILTERDATA}/${SP}_rev_merge_mod.bed

	for strtype in fwd rev
	do
	if [[ ${SAPNUM} == 3 ]];then
	cat ${FILTERDATA}/${SP}_${strtype}_R1.bed ${FILTERDATA}/${SP}_${strtype}_R2.bed ${FILTERDATA}/${SP}_${strtype}_R3.bed | cut -f 4 | sort -u | \
	awk -F"\t" 'NR==FNR{a[$1]=$1;next}{if(a[$4]){OFS="\t";print $0}}' - ${PCGGENE}.PCG.info.bed |
	bedtools intersect -nonamecheck -wo -s -a ${FILTERDATA}/${SP}_${strtype}_merge_mod.bed -b - | \
	awk -F"\t" '$13>=100{OFS="\t";print $1,$2,$3,$10,$11,$12}' - | sort -k1,1 -k2,2n - \
	> ${FILTERDATA}/${SP}_${strtype}_merge.bed
	elif [[ ${SAPNUM} == 2 ]];then
	cat ${FILTERDATA}/${SP}_${strtype}_R1.bed ${FILTERDATA}/${SP}_${strtype}_R2.bed | cut -f 4 | sort -u | \
	awk -F"\t" 'NR==FNR{a[$1]=$1;next}{if(a[$4]){OFS="\t";print $0}}' - ${PCGGENE}.PCG.info.bed |
	bedtools intersect -nonamecheck -wo -s -a ${FILTERDATA}/${SP}_${strtype}_merge_mod.bed -b - | \
	awk -F"\t" '$13>=100{OFS="\t";print $1,$2,$3,$10,$11,$12}' - | sort -k1,1 -k2,2n - \
	> ${FILTERDATA}/${SP}_${strtype}_merge.bed
	fi
	done

	for rep in `seq ${SAPNUM}`
	do
	for treat in IP Input
	do
	for strtype in fwd rev
	do
	samtools sort -@ 20 -o ${BAMPATH}/${SP}_${treat}${rep}_${strtype}_final_sort.bam ${BAMPATH}/${SP}_${treat}${rep}_${strtype}_final.bam
	if [[ ${SP} == "tae" ]] || [[ ${SP} == "tdi" ]] || [[ ${SP} == "ata" ]];then
	samtools index -c -@ 20 ${BAMPATH}/${SP}_${treat}${rep}_${strtype}_final_sort.bam
	cp ${BAMPATH}/${SP}_${treat}${rep}_${strtype}_final_sort.bam.csi ${BAMPATH}/${SP}_${treat}${rep}_${strtype}_final_sort.bam.bai
	else
	samtools index -@ 20 ${BAMPATH}/${SP}_${treat}${rep}_${strtype}_final_sort.bam
	fi
	done
	done
	done

	for strtype in fwd rev
	do
	cut -f 1-3 ${FILTERDATA}/${SP}_${strtype}_merge.bed > ${FILTERDATA}/${SP}_${strtype}_merge_cut.bed
	if [[ ${SAPNUM} == 3 ]];then
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${BAMPATH}/${SP}_IP1_${strtype}_final_sort.bam,Input1,${BAMPATH}/${SP}_Input1_${strtype}_final_sort.bam,${FILTERDATA}/${SP}_${strtype}_merge_cut.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${BAMPATH}/${SP}_IP2_${strtype}_final_sort.bam,Input2,${BAMPATH}/${SP}_Input2_${strtype}_final_sort.bam,${FILTERDATA}/${SP}_${strtype}_merge_cut.bed,bed
Sample3,Leaf,Leaf3,IP3,3,${BAMPATH}/${SP}_IP3_${strtype}_final_sort.bam,Input3,${BAMPATH}/${SP}_Input3_${strtype}_final_sort.bam,${FILTERDATA}/${SP}_${strtype}_merge_cut.bed,bed
" >${FILTERDATA}/DiffBind_${strtype}_conf.csv
	elif [[ ${SAPNUM} == 2 ]];then
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${BAMPATH}/${SP}_IP1_${strtype}_final_sort.bam,Input1,${BAMPATH}/${SP}_Input1_${strtype}_final_sort.bam,${FILTERDATA}/${SP}_${strtype}_merge_cut.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${BAMPATH}/${SP}_IP2_${strtype}_final_sort.bam,Input2,${BAMPATH}/${SP}_Input2_${strtype}_final_sort.bam,${FILTERDATA}/${SP}_${strtype}_merge_cut.bed,bed
" >${FILTERDATA}/DiffBind_${strtype}_conf.csv
	fi
	/usr/bin/Rscript scripts/diffbind2.R ${FILTERDATA}/DiffBind_${strtype}_conf.csv ${FILTERDATA}/${SP}_${strtype}_merge.bed ${FILTERDATA}/${SP}_${strtype}_merge_diffbind.txt
	done

	for strtype in fwd rev
	do
	if [[ ${SAPNUM} == 3 ]];then
	awk -F"\t" 'NR==FNR{a[$1"_"$2"_"$3]=$0;next}{print($0"\t"a[$1"_"$2"_"$3])}' \
	${FILTERDATA}/${SP}_${strtype}_merge_diffbind.txt ${FILTERDATA}/${SP}_${strtype}_merge.bed | \
	cut -f 1-6,10-15 | awk '{OFS="\t";$2+=1;print $0}' - >${FILTERDATA}/../${SP}_${strtype}.txt
	elif [[ ${SAPNUM} == 2 ]];then
	awk -F"\t" 'NR==FNR{a[$1"_"$2"_"$3]=$0;next}{print($0"\t"a[$1"_"$2"_"$3])}' \
	${FILTERDATA}/${SP}_${strtype}_merge_diffbind.txt ${FILTERDATA}/${SP}_${strtype}_merge.bed | \
	cut -f 1-6,10-13 | awk '{OFS="\t";$2+=1;print $0}' - >${FILTERDATA}/../${SP}_${strtype}.txt
	fi
	done

	mkdir -p data/${SP}/05Ratio
	cat ${FILTERDATA}/../${SP}_fwd.txt ${FILTERDATA}/../${SP}_rev.txt > data/${SP}/05Ratio/${SP}_gene_ratio.txt
	rm ${FILTERDATA}/../${SP}_fwd.bed ${FILTERDATA}/../${SP}_rev.bed ${FILTERDATA}/../${SP}.bed
	
	/usr/bin/Rscript scripts/samples_gene_merge.R ${GENEPATH}/raw Input ${GENEPATH}/${SP}_gene_expression.txt
	awk -F"\t" 'NR==FNR{a[$1]=$NF;next}{OFS="\t";if(a[$4]){print $0"\t"a[$4]}else{print $0"\t"0}}' \
	${GENEPATH}/${SP}_gene_expression.txt ${PCGGENE}.PCG.info.bed >${GENEPATH}/${SP}_gene_feature.txt
	
	echo "== ${SP} =="
	for rep in 1 2
	do
	samtools sort -@ 10 -o ${DATAPATH}/${SP}/03BED/${SP}_IP${rep}_final_sort.bam ${DATAPATH}/${SP}/03BED/${SP}_IP${rep}_final.bam
	cat ${FILTERDATA}/${SP}_fwd_R${rep}.bed ${FILTERDATA}/${SP}_rev_R${rep}.bed | \
	cut -f 1-3 - >${FILTERDATA}/../${SP}_R${rep}_cut.bed
	for treat in Input IP
	do
	if [[ ${SP} == "tae" ]] || [[ ${SP} == "tdi" ]] || [[ ${SP} == "ata" ]];then
	samtools index -c -@ 10 ${DATAPATH}/${SP}/03BED/${SP}_${treat}${rep}_final_sort.bam
	cp ${DATAPATH}/${SP}/03BED/${SP}_${treat}${rep}_final_sort.bam.csi ${DATAPATH}/${SP}/03BED/${SP}_${treat}${rep}_final_sort.bam.bai
	else
	samtools index -@ 10 ${DATAPATH}/${SP}/03BED/${SP}_${treat}${rep}_final_sort.bam
	fi
	done
	done
echo -e "SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller
Sample1,Leaf,Leaf1,IP1,1,${DATAPATH}/${SP}/03BED/${SP}_IP1_final_sort.bam,Input1,${DATAPATH}/${SP}/03BED/${SP}_Input1_final_sort.bam,${FILTERDATA}/../${SP}_R1_cut.bed,bed
Sample2,Leaf,Leaf2,IP2,2,${DATAPATH}/${SP}/03BED/${SP}_IP2_final_sort.bam,Input2,${DATAPATH}/${SP}/03BED/${SP}_Input2_final_sort.bam,${FILTERDATA}/../${SP}_R2_cut.bed,bed
Sample3,Leaf,Leaf3,IP3,3,${DATAPATH}/${SP}/03BED/${SP}_IP3_final_sort.bam,Input3,${DATAPATH}/${SP}/03BED/${SP}_Input3_final_sort.bam,${FILTERDATA}/../${SP}_R3_cut.bed,bed" >${FILTERDATA}/../DiffBind_corrplot.csv

done

