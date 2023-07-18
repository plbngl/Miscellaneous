#!/usr/bin/env Rscript

library('stringr')
library("data.table")

args = commandArgs(trailingOnly=TRUE)
id   = args[1]
file = paste0(id, ".somatic.p001.vep.tsv")

##### Read Varscan calls and filter for coverage, frequecny  and MAF #####
##########################################################################

ind     = read.table(file , header=T)  
ac      = str_split_fixed(ind$TUMOR_DP4 , "\\,",4) 
ac      = apply(ac, 2,as.numeric)
wc      = str_split_fixed(ind$NORMAL_DP4 , "\\,",4) 
wc      = apply(wc, 2,as.numeric)
ind$TUMOR_FREQ  = ind$TUMOR_AD / ind$TUMOR_DP
ind$NORMAL_FREQ = ind$NORMAL_AD / ind$NORMAL_DP

MIN_READS_ALT_FLAG_STRANDS_TUMOR  = ac[,3]>4 & ac[,4]>4
MIN_READS_ALT_FLAG_STRANDS_NORMAL = wc[,3]>4 & wc[,4]>4

frq_filt = ind$TUMOR_FREQ >0.333 & ind$NORMAL_FREQ <0.05
loh_filt = ind$TUMOR_GT == "1/1" & ind$TUMOR_FREQ >0.75 & ind$NORMAL_GT =="0/1" #& ind$NORMAL_FREQ < 0.66

ind$somatic_filters = frq_filt & MIN_READS_ALT_FLAG_STRANDS_TUMOR
ind$loh_filters     = loh_filt & MIN_READS_ALT_FLAG_STRANDS_TUMOR & MIN_READS_ALT_FLAG_STRANDS_NORMAL

n = max(str_count(ind$CSQ, "\\,"))+1

header = 'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|SIFT|PolyPhen|HGVS_OFFSET|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS'
header = unlist(strsplit(header, "|", fixed=T))

annot1  = sapply(ind$CSQ, function(x) strsplit(x, ",", fixed=T))
annot2  = lapply(annot1, function(x) data.frame(str_split_fixed(x, "\\|", length(header))))
test    = lapply(annot2, function(x) subset(x, x[,which(header=='PICK')]=='1'))
annot   = rbindlist(test)
colnames(annot)= header

ind2      = cbind(ind[,!(colnames(ind)%in% 'CSQ')],annot )
ind2      = as.data.frame(apply(ind2,2,function(x) gsub(" ", "", x)))

ind2$rare_absent_inPoP = as.numeric(ind2[,"MAX_AF"]) <0.005 | ind2[,"MAX_AF"] ==''

loh       = subset(ind2, ind2[,'loh_filters']==T)

ind2      = subset(ind2, rare_absent_inPoP==TRUE) ## remove SNPs for somatic / keep it for LOH

som = subset(ind2, ind2[,'somatic_filters']==T)

write.table(loh,   paste0( id, ".filtered.loh.annot.tsv"), sep="\t", quote=F, row.names=F)
write.table(som,   paste0( id, ".filtered.somatic.annot.tsv"), sep="\t", quote=F, row.names=F)


i = 0
for (df in list(som, loh)){
i=i+1
tp = c("somatic", "LOH")[i]  
    df = as.data.frame(df)

    
####### Plots (Manhattan style) ##########
##########################################
    
ch_file      = read.table("/home/paola/data/publicdata/genome/hg38_chrom_sizes")
ch_file[,1]  = paste0("chr", ch_file[,1])
chx          = cumsum(as.numeric(ch_file[,3]))
delta        = chx - ch_file[,3]
names(delta) = ch_file[,1]
midpoint     = chx - ch_file[,3]/2
    
    
colors        = c('red', "orange", "yellow", 'green')
names(colors) = c('HIGH',  'MODERATE',   'LOW',  'MODIFIER' )
mah           = data.frame(df[,c('SYMBOL', 'CHROM', 'POS', 'IMPACT', 'VARIANT_CLASS', 'rare_absent_inPoP')])
mah$chr       = gsub("chr", '', mah$CHROM)
mah$chr       = gsub("X", '23', mah$chr)
mah$chr       = gsub("Y", '24', mah$chr)
mah$freq_mut  = as.numeric(df$TUMOR_AD)/(as.numeric(df$TUMOR_AD) + as.numeric(df$TUMOR_RD))
mah           = mah[order( as.numeric(mah$chr), mah$POS ),]
mah$offset    = delta[mah$CHR]
mah$genpos    = as.numeric(mah$POS) + mah$offset
    
    
pdf( paste0( id, ".",tp,".Manhattan_plot_varscan.pdf"), height = 4, width = 8)
par(mar=c(4, 2, 4, 6)) 
plot(mah$genpos, mah$freq_mut, xlim=c(0,max(chx)),axes=F , pch= c(24,21)[(mah$VARIANT_CLASS=="SNV")+1],
     ylim=c(0.3,1),
     main=id,
     bg=colors[mah$IMPACT], col = c("black", "white") [mah$rare_absent_inPoP+1],
     ylab = "Allele.freq", xlab=NA)
axis(side = 2, at=c(0,0.25,0.5,0.75,1),lwd = 2, line=-1)
abline(v=chx, col="gray", lty=2)
abline(h=0.5, col="blue")
axis(side = 1, lwd=2,at =c(0, midpoint, max(chx)),
     labels =c(NA,ch_file[,1] ,NA), las=2)


lbls = subset(mah,IMPACT!='MODIFIER' & rare_absent_inPoP==TRUE)
if(nrow(lbls)>0){
  text(lbls$genpos, lbls$freq_mut, labels = lbls$SYMBOL, pos = 4, srt=45,xpd=T)
}

    if(nrow(mah)>1){
ag = aggregate( genpos~ IMPACT,mah, length)
legend("topright", inset=c(-0.15,-0.15), col= colors[ag[,1]],cex=0.8,bty = 'n',
       legend=paste0(names(colors[ag[,1]]),"\n(", ag[,2], ")\n"), pch=15, title="Impact",xpd=TRUE)
}else{
        
  legend("topright", inset=c(-0.15,-0.15), col= colors[mah[,"IMPACT"]],cex=0.8,bty = 'n',
       legend=paste0(names(colors[mah[,"IMPACT"]]),"\n(", 1, ")\n"), pch=15, title="Impact",xpd=TRUE)
      
    }

legend("bottomright", inset=c(-0.08,0),cex=0.8,bty = 'n',
       legend=c("snv", "indel"), pch=c(16,17),xpd=TRUE)
dev.off()


colpal = c(rainbow(10)[c(1,2,4, 7,9,10)], 'yellow', 'gray' , 
                             cm.colors(2), terrain.colors(4),1:6)


if(nrow(df)>0){

tb=table(df$Consequence)
#tb=sort(tb,decreasing = T)
names(tb) = paste0(names(tb)," (" ,tb, ")")


pdf( paste0( id,".", tp,".Piechart_varscan.pdf"), height = 4, width = 9)
par (mfrow=c(1,2),xpd=TRUE, mar=c(4,0,4,0))
pie(tb, col=colpal,border="white", labels = NA, radius=1.1, main = id)
plot.new()
legend("topright", fill =colpal,legend = names(tb) ,border = NA, bty="n")
dev.off()
    }
}
 
