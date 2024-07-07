mkdir /project2/xinhe/1kg/1000G_Phase3_plink/b38
cd /project2/xinhe/1kg/1000G_Phase3_plink/b38

# Download raw files (b38 from https://www.cog-genomics.org/plink/2.0/resources)
# The links in here may be changed in future
# "-O" to specify output file name
wget -O all_phase3.pgen.zst "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1"
wget -O all_phase3.pvar.zst "https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz&dl=1"
wget -O all_phase3.psam "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1"

# Use PLINK2 to decompress the pgen file
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

# "vzs" modifier to directly operate with pvar.zst
# "--chr 1-22" select variants on chr1:22
# "--output-chr 26" uses numeric chromosome codes
# "--maf 0.01" limits to MAF > 0.01
# "--max-alleles 2": filters out the multiallelic variants,
# PLINK 1 binary does not allow multiallelic variants
# "--rm-dup" removes duplicate-ID variants
# "--set-missing-var-id" replaces missing IDs with a pattern
# restricting to autosomal SNPs with MAF>0.01
plink2 --pfile all_phase3 vzs \
--chr 1-22 \
--allow-extra-chr \
--output-chr 26 \
--maf 0.01 \
--max-alleles 2 \
--rm-dup exclude-all \
--set-missing-var-ids '@:#_\$a_\$r' \
--make-pgen \
--out all_phase3_autosomes

# Prepare sub-population filter file
awk 'NR == 1 || $5 == "EUR" {print $1}' all_phase3.psam > EUR_1kg_samples.txt

awk 'NR == 1 || $5 == "EAS" {print $1}' all_phase3.psam > EAS_1kg_samples.txt

awk 'NR == 1 || $5 == "AFR" {print $1}' all_phase3.psam > AFR_1kg_samples.txt

# Generate sub-population fileset
# Convert the 1000 Genomes data to PLINK 1 binary format
plink2 --pfile all_phase3_autosomes \
--keep EUR_1kg_samples.txt \
--make-bed \
--out 1000G_EUR_phase3_autosomes

plink2 --pfile all_phase3_autosomes \
--keep EAS_1kg_samples.txt \
--make-bed \
--out 1000G_EAS_phase3_autosomes

plink2 --pfile all_phase3_autosomes \
--keep AFR_1kg_samples.txt \
--make-bed \
--out 1000G_AFR_phase3_autosomes

# Split bed/bim/fam by chromosome
mkdir 1000G_EUR_phase3_b38
for i in {1..22}
do
plink2 --bfile 1000G_EUR_phase3_autosomes \
--chr $i \
--make-bed \
--out 1000G_EUR_phase3_b38/1000G_EUR_phase3_b38_$i
done

mkdir 1000G_EAS_phase3_b38
for i in {1..22}
do
plink2 --bfile 1000G_EAS_phase3_autosomes \
--chr $i \
--make-bed \
--out 1000G_EAS_phase3_b38/1000G_EAS_phase3_b38_$i
done

mkdir 1000G_AFR_phase3_b38
for i in {1..22}
do
plink2 --bfile 1000G_AFR_phase3_autosomes \
--chr $i \
--make-bed \
--out 1000G_AFR_phase3_b38/1000G_AFR_phase3_b38_$i
done
