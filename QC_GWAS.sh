plink --bfile ${1} --missing 
plink --bfile ${1} --geno 0.1 --make-bed --out ${1}_2 
plink --bfile ${1}_2 --mind 0.1 --make-bed --out ${1}_3 
plink --bfile ${1}_3 --geno 0.02 --make-bed --out ${1}_4 
plink --bfile ${1}_4 --mind 0.02 --make-bed --out ${1}_5 

plink --bfile ${1}_5 --freq --out MAF_check 
plink --bfile ${1}_5 --maf 0.05 --make-bed --out ${1}_8 
plink --bfile ${1}_8 --hwe 1e-6 --hwe-all --make-bed --out ${1}_9 

#final_step





