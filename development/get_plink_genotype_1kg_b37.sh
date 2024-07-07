mkdir /project2/xinhe/1kg/1000G_Phase3_plink/b37
cd /project2/xinhe/1kg/1000G_Phase3_plink/b37

# Download raw files (b37 from https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3)
# The links in here may be changed in future
# "-O" to specify output file name
wget -O all_phase3.pgen.zst "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
wget -O all_phase3.pvar.zst "https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1"
wget -O all_phase3.psam "https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?rlkey=0yyifzj2fb863ddbmsv4jkeq6&dl=1"

# Use PLINK2 to decompress the pgen and pvar files
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

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

# optionally:
# Download genetic maps and insert them using PLINK1.9
wget https://www.dropbox.com/s/slchsd0uyd4hii8/genetic_map_b37.zip
unzip genetic_map_b37.zip

plink --bfile 1000G_EUR_phase3_autosomes \
--cm-map genetic_map_b37/genetic_map_chr@_combined_b37.txt \
--make-bed \
--out 1000G_EUR_phase3_autosomes_b37

plink --bfile 1000G_EAS_phase3_autosomes \
--cm-map genetic_map_b37/genetic_map_chr@_combined_b37.txt \
--make-bed \
--out 1000G_EAS_phase3_autosomes_b37

plink --bfile 1000G_AFR_phase3_autosomes \
--cm-map genetic_map_b37/genetic_map_chr@_combined_b37.txt \
--make-bed \
--out 1000G_AFR_phase3_autosomes_b37

# Split bed/bim/fam by chromosome
mkdir 1000G_EUR_phase3_b37
for i in {1..22}
do
plink2 --bfile 1000G_EUR_phase3_autosomes_b37 \
--chr $i \
--make-bed \
--out 1000G_EUR_phase3_b37/1000G_EUR_phase3_b37_chr$i
done

mkdir 1000G_EAS_phase3_b37
for i in {1..22}
do
plink2 --bfile 1000G_EAS_phase3_autosomes_b37 \
--chr $i \
--make-bed \
--out 1000G_EAS_phase3_b37/1000G_EAS_phase3_b37_chr$i
done

mkdir 1000G_AFR_phase3_b37
for i in {1..22}
do
plink2 --bfile 1000G_AFR_phase3_autosomes_b37 \
--chr $i \
--make-bed \
--out 1000G_AFR_phase3_b37/1000G_AFR_phase3_b37_chr$i
done

