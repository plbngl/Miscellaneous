library(DNAcopy)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
sample_name   = args[1]
cnvkit_dir    = args[2]
variants_dir  = args[3]

setwd(cnvkit_dir)

ch_file = read.table("/home/paola/references/hg38_chrom_sizes")
ch_file[,1]=paste0("chr", ch_file[,1])
chx = cumsum(as.numeric(ch_file[,3]))
delta = chx - ch_file[,3]
names(delta) = ch_file[,1]
midpoint = chx - ch_file[,3]/2
chrom_numbers = setNames( 1:24, paste0("chr", c(1:22, "X", "Y")))

cn <- read.table(sprintf("%s.cnr",sample_name ), header=T)
cn = subset(cn, depth>0)
cs <- read.table(sprintf("%s.cns",sample_name ), header=T)
CNA.object <-CNA( genomdat = cn[,'log2'], 
                   chrom = cn[,'chromosome'],
                 maploc = cn[,'start'], data.type = 'logratio')

CNA.smoothed <- smooth.CNA(CNA.object)
CNA.smoothed$pos = as.numeric(CNA.smoothed$maploc) + delta[CNA.smoothed$chrom]
mah = cs
mah$chr = chrom_numbers[mah$chromosome]
mah = mah[order( as.numeric(mah$chr), mah$start ),]
mah$offset = delta[mah$chromosome]
mah$genpos = as.numeric(mah$start) + mah$offset
mah$genend = as.numeric(mah$end) + mah$offset

colvect=c("gray", "gray20")
colvect2=c("orange", "red")

######    Plot 1 Plot separate chromosome - no scale    ################
########################################################################

pdf(sprintf("%s_CNVkit_byCHR.pdf", sample_name), height = 10, width = 15)
par(mfrow=c(6,4), mar=c(1,4,1,1))
for (chromo in 1:24){
lfc = CNA.smoothed[CNA.smoothed$chrom== names(which(chrom_numbers==chromo)),]
segm = mah[mah$chr==chromo,]
plot(x=lfc$maploc,y=lfc$Sample.1, pch=".", 
     col = colvect  [2],
     axes=F, ylim=c(-2.5,2.5), 
     xlab="", ylab="",
    # xlim=c(1,ch_file[1,3])
    )
axis(side = 2,lwd = 2, line=0) 
axis(side = 1, lwd=2,at =c(0, ch_file[chromo,3]))
mtext(ch_file[chromo,1], line=-1)
mtext("Copy ratio (log2)", side=2, line=2, cex=0.7)
segments(segm$start, segm$log2, x1 = segm$end, 
         col = colvect2[(abs(segm$log2)>0.5) +1],  lwd = 2)
    }
dev.off()

##### Add B allele frequency #######
####################################

vars = read.table(paste0(variants_dir, sample_name, ".exome_snps.tsv"), header=F) 
vars = subset(vars, V6>20)
tabs = str_split_fixed(vars[,7], "\\,", 4)
tabs = apply(tabs, 2, as.numeric)
VAF = rowSums(tabs[,3:4])/ rowSums(tabs)
df1 = data.frame( chr =vars[,1] , pos = as.numeric(vars$V2), vaf = as.numeric(VAF), id =paste(vars[,1], vars[,2], sep=":"))

extras = read.table(paste0(variants_dir, sample_name, ".germline.variants.tsv"), header=F) 
extras = subset(extras, V6>20)
tabs   = str_split_fixed(extras[,8], ",|=", 5)
tabs   = apply(tabs[,2:5], 2, as.numeric)
VAF2   = rowSums(tabs[,3:4])/ rowSums(tabs)
df2    = data.frame(chr =extras[,1] , pos = as.numeric(extras$V2), vaf = as.numeric(VAF2), id =paste(extras[,1], extras[,2], sep=":"))


df = rbind(df1,df2)
df = df[!duplicated(df$id),]
df$pos_abs = df$pos + delta[df$chr]

colvect3 = c("purple", "blue")
colvect2 = c("blue", "red")

### Plot 2 Plot separate chromosome - no scale and B allele frequecny 
######################################################################## 

pdf(sprintf("%s_CNVkit_byCHRand_denseBAF.pdf", sample_name), height = 10, width = 18)
par(mfrow=c(6,4), mar=c(1,3,1,3))
for (chromo in 1:24){ 
lfc = CNA.smoothed[CNA.smoothed$chrom== names(which(chrom_numbers==chromo)),]
segm = mah[mah$chr==chromo,]
plot(x=lfc$maploc,y=lfc$Sample.1, pch=".", 
     col = 'red',
     axes=F, ylim=c(-2.5,2.5), 
     xlab="", ylab="",
    # xlim=c(1,ch_file[1,3])
    )
axis(side = 2,lwd = 2, line=-0.8) 
mtext(ch_file[chromo,1], line=0)
mtext("Copy ratio (log2)", side=2, line=1.5, cex=0.75)

segments(segm$start, segm$log2, x1 = segm$end, 
         col = colvect2[(abs(segm$log2)>0.5) +1],  lwd = 2)
par(new=TRUE)

dfc = subset(df, chr == names(which(chrom_numbers==chromo)))

plot(dfc[,2:3], pch=".", 
     col = 'blue',
     axes=F, ylim=c(0,1), 
     xlab="", ylab="",
    # xlim=c(1,ch_file[1,3])
    )
axis(side = 4,lwd = 2, line=-0.8 )
mtext("VAF", side=4, line=1.5, cex=0.75)
# par(new=TRUE)
# plot(x=lfc$maploc,y=lfc$Sample.1, type="n",  
#      axes=F, ylim=c(-2.5,2.5), 
#      xlab="", ylab="",
#     )
# segments(segm$start, segm$log2, x1 = segm$end, 
#          col = colvect2[(abs(segm$log2)>0.75) +1],  lwd = 2)

}
dev.off()


### just use the panel for the following plots
df = df1
df$pos_abs = df$pos + delta[df$chr]



### Plot 3 Plot separate chromosome - scaled and B allele frequecny 
########################################################################

pdf(sprintf("%s_CNVkit_karyoview.pdf",sample_name), height = 18, width = 9)
par(mfrow=c(24,1), mar = c(1,5,0,1))
for (chromo in 1:24){ 
lfc = CNA.smoothed[CNA.smoothed$chrom== names(which(chrom_numbers==chromo)),]
segm = mah[mah$chr==chromo,]
plot(x=lfc$maploc,y=lfc$Sample.1, pch=".", 
     col ='red',
     axes=F, ylim=c(-2.5,2.5), xlab="", ylab="",
     xlim=c(1,ch_file[1,3]))
axis(side = 2,lwd = 1.5, line=-2)
axis(side = 1, lwd=1.5,at =c(0, ch_file[chromo,3]))
mtext(ch_file[chromo,1],side=2, line=1, las=2)
mtext("Copy ratio", side=2, line=0, cex=0.7)

segments(segm$start, segm$log2, x1 = segm$end, 
         col = colvect2[(abs(segm$log2)>0.5) +1],  lwd = 2)

par(new=TRUE)

dfc = subset(df, chr == names(which(chrom_numbers==chromo)))

plot(dfc[,2:3], pch=".", 
     col = 'blue',
     axes=F, ylim=c(0,1), 
     xlab="", ylab="",
       xlim=c(1,ch_file[1,3])
    )

}
dev.off()

### Plot 4 Whole Genome view cnvkit and B allele frequecny 
########################################################################
#options(repr.plot.width=15, repr.plot.height=8)
pdf(sprintf("%s_WG_plot.pdf", sample_name), height=5,width = 10)
par(mfrow=c(2,1), mar = c(4,4,1,1))
plot(x=CNA.smoothed$pos,y=CNA.smoothed$Sample.1, pch=".",  
     col = colvect[(chrom_numbers[CNA.smoothed$chrom] %% 2)+1],
     axes=F, ylim=c(-3,3), xlab="", ylab="Copy ratio (log2)")
axis(side = 2,lwd = 2, line=-1, las=2)
axis(side = 1, lwd=2,at =c(0, midpoint, max(chx)),srt=45,
     labels =c(NA,ch_file[,1] ,NA), las=2)

segments(mah$genpos, mah$log2, x1 = mah$genend, 
         col = c('orange', 'red')[(abs(mah$log2)>0.75) +1],  lwd = 2)
mtext(sample_name, line=-1, font=2)


plot(df$pos_abs,df$vaf, pch=".", 
      col = colvect3  [(chrom_numbers[df$chr] %% 2)+1],
     axes=F, ylim=c(0,1), 
     xlab="", ylab="B allele frequency"
    )
axis(side = 2,lwd = 2, line=-1, las=2) 
axis(side = 1, lwd=2,at =c(0, midpoint, max(chx)),srt=45,
     labels =c(NA,ch_file[,1] ,NA), las=2)
dev.off()

