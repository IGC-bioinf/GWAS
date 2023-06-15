## Right calling first(WARNING PARALLEL 40)
## bash script.sh /path/to/ref /path/to/bed

for i in $(ls ./reads/*.fastq.gz | sed "s/_R1.fastq.gz//g" | sed "s/_R2.fastq.gz//g" | sort | uniq)
do
	bwa mem -t 40 ref_37.fna ${i}_R1.fastq.gz ${i}_R2.fastq.gz | samtools view -@ 40 -Sb  > ${i}_without_RG.bam
	java -jar gatk-package-4.1.1.0-local.jar AddOrReplaceReadGroups  -I ${i}_without_RG.bam -O ${i}.bam -SORT_ORDER coordinate  --RGPL ILLUMINA --RGSM ${i} --RGPU exome_seq --RGLB lib_${i}
	samtools index ${i}.bam
	rm ${i}_without_RG.bam
done
task(){
	java -jar gatk-package-4.1.1.0-local.jar HaplotypeCaller --native-pair-hmm-threads 1  -ERC GVCF -L /PATH/TO/BED  -R ref_37.fna -I $1 -O ${1}.g.vcf
}

N=40
for thing in $(ls *.bam | grep -v "RG" | grep -v "Undetermined.bam")
do
	((i=i%N)); ((i++==0)) && wait
	task "$thing" & 
done
java -jar gatk-package-4.1.1.0-local.jar CombineGVCFs $(find ./ -name *g.vcf | sed "s~"^./"~ \ -V ./~g") -R ref_37.fna -O Combined_non_geno.g.vcf
java -jar gatk-package-4.1.1.0-local.jar GenotypeGVCFs \
   -R ref_37.fna \
   -V Combined_non_geno.g.vcf \
   -O GENOTYPED_no_ann.vcf.gz
bcftools annotate --threads 8 -x ID -I +'%CHROM:%POS:%REF:%ALT' -Oz -o GENOTYPED.vcf.gz GENOTYPED_no_ann.vcf.gz


#files need to be in directory: splitted 1000G, target vcf.gz file. 
#For GWAS: sex.txt, family.txt, parents.txt, pheno.txt.


mkdir output_isec
for i in {1..22}
do
	bcftools isec --threads 40 -p ./output_isec/chr${i} GENOTYPED.vcf.gz ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz -Oz
done
bcftools isec --threads 40 -p ./output_isec/chrX GENOTYPED.vcf.gz ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz -Oz
bcftools isec --threads 40 -p ./output_isec/chrY GENOTYPED.vcf.gz ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz -Oz
bcftools concat -o 1000G_isec_target_no_ann.vcf $(find ./ -name 0003.vcf.gz)
bcftools annotate --threads 8 -x ID -I +'%CHROM:%POS:%REF:%ALT' -Oz -o 1000G_isec_target.vcf.gz 1000G_isec_target_no_ann.vcf
####STARTING GWAS#####


# Covert to bed
plink --vcf $1 --make-bed --out proj_1
# Update sex, family, pheno

plink --bfile proj_1 --update-sex sex.txt --make-bed --out proj_1
plink --bfile proj_1 --update-ids family.txt --make-bed --out proj_1
plink --bfile proj_1 --update-parents parents.txt --make-bed --out proj_1
plink --bfile proj_1 --pheno pheno.txt --make-bed --out proj_1

#QC

plink --bfile proj_1 --missing
Rscript --no-save hist_miss.R
plink --bfile proj_1 --geno 0.2 --make-bed --out proj_2
plink --bfile proj_2 --mind 0.2 --make-bed --out proj_3
plink --bfile proj_3 --geno 0.02 --make-bed --out proj_4
plink --bfile proj_4 --mind 0.02 --make-bed --out proj_5
Rscript --no-save gender_check.R
grep "PROBLEM" plink.sexcheck| awk '{print$1,$2}'> sex_discrepancy.txt
plink --bfile proj_5 --remove sex_discrepancy.txt --make-bed --out proj_6
plink --bfile proj_5 --impute-sex --make-bed --out proj_6
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' proj_6.bim > snp_1_22.txt
plink --bfile proj_6 --extract snp_1_22.txt --make-bed --out proj_7
plink --bfile proj_7 --freq --out MAF_check
plink --bfile proj_7 --maf 0.05 --make-bed --out proj_8
plink --bfile proj_8 --hwe 1e-10 --hwe-all --make-bed --out proj_9
#plink --bfile proj_9 --exclude inversion.txt --range --indep-pairwise 50 5 0.2 --out indepSNP
plink --bfile proj_9 --extract indepSNP.prune.in --het --out R_check
Rscript --no-save check_heterozygosity_rate.R
Rscript --no-save heterozygosity_outliers_list.R
sed 's/"// g' fail-het-qc.txt | awk '{print$1, $2}'> het_fail_ind.txt
plink --bfile proj_9 --remove het_fail_ind.txt --make-bed --out proj_10
plink --bfile proj_10 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
Rscript --no-save Relatedness.R
plink --bfile proj_10 --filter-founders --make-bed --out proj_11
plink --bfile proj_11 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
plink --bfile proj_11 --missing



#Stratification
plink --vcf 1000G_isec_target.vcf.gz --make-bed --out ALL.2of4intersection.20100804.genotypes
plink --bfile ALL.2of4intersection.20100804.genotypes --geno 0.2 --allow-no-sex --make-bed --out 1kG_MDS
plink --bfile 1kG_MDS --mind 0.2 --allow-no-sex --make-bed --out 1kG_MDS2
plink --bfile 1kG_MDS2 --geno 0.02 --allow-no-sex --make-bed --out 1kG_MDS3
plink --bfile 1kG_MDS3 --mind 0.02 --allow-no-sex --make-bed --out 1kG_MDS4
plink --bfile 1kG_MDS4 --maf 0.05 --allow-no-sex --make-bed --out 1kG_MDS5
awk '{print$2}'  proj_11.bim > HapMap_SNPs.txt
plink --bfile 1kG_MDS5 --extract HapMap_SNPs.txt --make-bed --out 1kG_MDS6
awk '{print$2}' 1kG_MDS6.bim > 1kG_MDS6_SNPs.txt
plink --bfile proj_11 --extract 1kG_MDS6_SNPs.txt --recode --make-bed --out HapMap_MDS
awk '{print$2,$4}' HapMap_MDS.map > buildhapmap.txt
plink --bfile 1kG_MDS6 --update-map buildhapmap.txt --make-bed --out 1kG_MDS7
awk '{print$2,$5}' 1kG_MDS7.bim > 1kg_ref-list.txt
plink --bfile HapMap_MDS --reference-allele 1kg_ref-list.txt --make-bed --out HapMap-adj
awk '{print$2,$5,$6}' 1kG_MDS7.bim > 1kGMDS7_tmp
awk '{print$2,$5,$6}' HapMap-adj.bim > HapMap-adj_tmp
sort 1kGMDS7_tmp HapMap-adj_tmp |uniq -u > all_differences.txt
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt
plink --bfile HapMap-adj --flip flip_list.txt --reference-allele 1kg_ref-list.txt --make-bed --out corrected_hapmap
awk '{print$2,$5,$6}' corrected_hapmap.bim > corrected_hapmap_tmp
sort 1kGMDS7_tmp corrected_hapmap_tmp |uniq -u  > uncorresponding_SNPs.txt
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
plink --bfile corrected_hapmap --exclude SNPs_for_exlusion.txt --make-bed --out HapMap_MDS2
plink --bfile 1kG_MDS7 --exclude SNPs_for_exlusion.txt --make-bed --out 1kG_MDS8
plink --bfile HapMap_MDS2 --bmerge 1kG_MDS8.bed 1kG_MDS8.bim 1kG_MDS8.fam --allow-no-sex --make-bed --out MDS_merge2
plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2
plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2
awk '{print$1,$1,$2}' integrated_call_samples_v3.20130502.ALL.panel > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt
awk '{print$1,$2,"BEL"}' HapMap_MDS.fam>racefile_own.txt
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt
Rscript MDS_merged.R
awk '{ print $1,$2 }' MDS_merge2.mds > EUR_MDS_merge2
plink --bfile proj_11 --allow-no-sex --keep EUR_MDS_merge2 --make-bed --out proj_13
plink --bfile proj_13 --extract indepSNP.prune.in --genome --out proj_13
plink --bfile proj_13 --read-genome proj_13.genome --cluster --mds-plot 10 --out proj_13_mds
awk '{print$1, $2, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' proj_13_mds.mds > covar_mds.txt


#final_step
plink --bfile proj_13 --assoc --out assoc_results
plink --bfile proj_13 --covar covar_mds.txt --logistic --hide-covar --out logistic_results
awk '!/'NA'/' logistic_results.assoc.logistic > logistic_results.assoc_2.logistic
plink --bfile proj_13 -assoc --adjust --out adjusted_assoc_results
awk '{ if ($4 >= 11595000 && $4 <= 21605000) print $2 }' proj_13.bim > subset_snp_chr_22.txt
plink --bfile proj_13 --extract subset_snp_chr_22.txt --make-bed --out HapMap_subset_for_perm
head sorted_subset.txt

Rscript --no-save Manhattan_plot.R
Rscript --no-save QQ_plot.R

## For PRSrice ask roma
