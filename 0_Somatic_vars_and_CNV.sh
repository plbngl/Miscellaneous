!/bin/bash


samp=$1
dir=$2
WTBAM=$3

refGenome="/home/paola/references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
VARSCAN='/home/paola/software/VarScan.v2.3.9.jar'
VEP='/home/paola/software/ensembl-vep/vep'


mkdir -p $samp
cd $samp

########## Call Somatic Variants using Varscan ############
###########################################################

if [ "$4" != "--novarscan" ]
then

        pileup="samtools mpileup -f $refGenome -q 1 -B -r {} \
        $WTBAM $dir/$samp/enrich_seq/$samp.bam | java -jar $VARSCAN somatic \
        --mpileup $samp.somatic.{} --output-vcf"


	echo "Starting Varscan calling sample $samp"

	cat /home/paola/references/hg38_chrom_sizes | awk 'FS="\t", OFS="" {print "chr"$1,":"$2"-"$3}' | parallel -j 48 $pileup

	for CH in {1..22} X Y; do echo $samp.somatic.chr$CH:*.indel.vcf; done > list.txt
	vcf-concat -f list.txt > $samp.somatic.concat.indel.vcf

	for CH in {1..22} X Y; do echo $samp.somatic.chr$CH:*.snp.vcf; done > list.txt
	vcf-concat -f list.txt > $samp.somatic.concat.snp.vcf

	rm $samp.somatic.chr*

fi

#### NOTE ####
## This by chromosome approach does not consider the extra contigs (chrUN, HLA...)- which are ~50% of all variants if called on the whole genome

bcftools view $samp.somatic.concat.indel.vcf --include 'INFO/SPV < 0.001' > $samp.somatic.concat.indel.p001.vcf
bcftools view $samp.somatic.concat.snp.vcf --include 'INFO/SPV < 0.001' > $samp.somatic.concat.snp.p001.vcf
vcf-concat $samp.somatic.concat.snp.p001.vcf $samp.somatic.concat.indel.p001.vcf > $samp.somatic.concat.all.p001.vcf

echo "Starting VEP annotation"

$VEP --cache -i $samp.somatic.concat.all.p001.vcf \
-o $samp.somatic.p001.vep.vcf \
--vcf --max_af --regulatory --sift b --polyphen b \
--variant_class --numbers --symbol \
--canonical --biotype --tsl --check_existing \
--fork 4 --force_overwrite --pick_order biotype,ccds,rank,canonical \
--flag_pick --fasta $refGenome --hgvs \
--custom /home/paola/references/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN


vk vcf2tsv wide --print-header $samp.somatic.p001.vep.vcf > $samp.somatic.p001.vep.tsv

Rscript --vanilla 1_filter_and_plot.R $samp

echo "Finish Variant Calling."
echo "Star CNV Calling."


##### CNV CALLING CNVKIT #####
##############################

BAM=$dir/$samp/enrich_seq/$samp.bam 

#cnvkit batch $BAM --normal $WTBAM \
#    --targets /home/paola/references/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.no.anno \
#    --fasta $refGenome  --output-reference reference.cnn \
#    --output-dir cnvkit_results/ \
#    --scatter --drop-low-coverage -p 48 

cnvkit batch $BAM \
    -r ../reference.cnn \
    --output-dir cnvkit_results/ \
    --scatter --drop-low-coverage -p 48 


##### Calculate B allele frequency #####
########################################

##### Call common exonic variants (dbSNP > 1%) including homo ref

pileup="bcftools mpileup -Ou -f $refGenome -q 1 -B -R /home/paola/references/exome_snp_positions/chr{}.exome.snp.pos \
         $dir/$samp/enrich_seq/$samp.bam | bcftools call -m  -Ou | \
        bcftools query -i 'QUAL>20 && DP>8 && AN==2' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/DP\t%INFO/DP4\n' > $samp.exome_snps.chr{}.tsv" 


parallel -j 48 $pileup ::: {1..22} X Y

cat $samp.exome_snps.chr* > $samp.exome_snps.tsv
rm $samp.exome_snps.chr*

##### Select also the germline variants (only non referece) previously called by varscan

bcftools query -i 'GPV < 0.01 && SPV > 0.05' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/DP[\t%SAMPLE=%DP4]\n' $samp.somatic.concat.snp.vcf > germline.snps.tsv 
bcftools query -i 'GPV < 0.01 && SPV > 0.05' -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/DP[\t%SAMPLE=%DP4]\n' $samp.somatic.concat.indel.vcf > germline.indels.tsv 
cat germline.snps.tsv  germline.indels.tsv  >  $samp.germline.variants.tsv
rm germline.snps.tsv
rm germline.indels.tsv

##### Plot results #######

Rscript 2_CNV_plots.R  $samp $PWD/cnvkit_results/ $PWD/
cd ..
mkdir -p cnv_summaries
ln -s $PWD/$samp/cnvkit_results/${samp}_WG_plot.pdf cnv_summaries/
ln -s $PWD/$samp/cnvkit_results/${samp}_CNVkit_karyoview.pdf cnv_summaries/


echo "Finish CNV Calling."

