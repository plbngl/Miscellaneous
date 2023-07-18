suppressPackageStartupMessages(library(SoupX))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(DoubletFinder))
suppressPackageStartupMessages(library(repr))
suppressPackageStartupMessages(library(ggpubr))
library(scDblFinder)

args = commandArgs(trailingOnly=TRUE)
sampleid <- args[1]

######## Read individual CellRanger Outputs,     ########
######## perform hard filters and normalization, ########
######## doublet detectio nand soupX correction  ########

rawdir = sprintf('/mnt/vol-hdd/projects/ink_scrna/cellranger_output/%s/outs/multi/count/', sampleid)
outdir = sprintf('/mnt/vol-hdd/projects/ink_scrna/seurat_output/preprocessing/%s/', sampleid)

dir.create(outdir )
print(paste("Start processing:",Sys.time()))
fdata =Read10X(file.path(rawdir, "raw_feature_bc_matrix"))
adata_full <- CreateSeuratObject(counts = fdata, project =sampleid )
adata_full[["percent.mt"]] <- PercentageFeatureSet(adata_full, pattern = "^MT-")

### Fist doublet detection with scDblFinder to be performed on lightly filtered data

adata1  <- subset( adata_full,
    nFeature_RNA > 200 &  
    nCount_RNA >1000
    ) 


sce <- scDblFinder(GetAssayData(object = adata1, slot = "counts"))
sce_results = data.frame(SummarizedExperiment::colData(sce))
adata1@meta.data = cbind(adata1@meta.data,sce_results)

######################################
### Apply standard hard filters and QC plots
######################################
adata <- subset( adata1,
    nFeature_RNA > 200 &  
    nFeature_RNA < 6000 & ## potential doublets hard filter
    percent.mt < 10 &
    nCount_RNA >1000
    ) 

#### Additional useful QC measures
adata[rownames(adata) != "MALAT1",] -> adata.nomalat
apply( adata.nomalat@assays$RNA@counts,  2,  max) -> adata.nomalat$largest_count
apply( adata.nomalat@assays$RNA@counts, 2, which.max) -> adata.nomalat$largest_index
rownames(adata.nomalat)[adata.nomalat$largest_index] -> adata.nomalat$largest_gene
100 * adata.nomalat$largest_count / adata.nomalat$nCount_RNA -> adata.nomalat$percent.Largest.Gene
adata.nomalat$largest_gene -> adata$largest_gene
adata.nomalat$percent.Largest.Gene -> adata$percent.Largest.Gene
rm(adata.nomalat)

options(repr.plot.width=12, repr.plot.height=7)
vp <-VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 'percent.Largest.Gene'), ncol = 4)

plot1 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mt",)+ NoLegend()
plot2 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ NoLegend()
plot3 <-FeatureScatter(adata,feature1 = "nCount_RNA", feature2 = "percent.Largest.Gene")+ NoLegend()

png(paste0(outdir, sampleid, ".postfilterQC.png"), width = 700, height = 576)
ggarrange(vp, plot1 + plot2 + plot3, ncol = 1, nrow = 2, common.legend = TRUE,legend="none")
dev.off()

###################
### Some other metrics for complexity and high expr genes ## optional
###################
as_tibble(adata[[]], rownames="Cell.Barcode") -> qc.metrics
qc.metrics %>% mutate(complexity=log10(nFeature_RNA) / log10(nCount_RNA))  -> qc.metrics
lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA)) -> complexity.lm
qc.metrics %>% mutate(
 complexity_diff = log10(nFeature_RNA) - ((log10(qc.metrics$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
  ) -> qc.metrics
p1<- qc.metrics %>%
  ggplot(aes(x=complexity_diff)) +
  geom_density(fill="green")

min(c(max(qc.metrics$complexity_diff),0-min(qc.metrics$complexity_diff))) -> complexity_scale

p2<- qc.metrics %>%
  mutate(complexity_diff=replace(complexity_diff,complexity_diff< -0.1,-0.1)) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) +
  geom_point(size=0.5) +
  geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) +
  scale_colour_gradient2(low="blue2",mid="grey",high="red2")
qc.metrics %>%
  group_by(largest_gene) %>%
  count() %>%
  arrange(desc(n)) -> largest_gene_list
largest_gene_list %>%
  filter(n>140) %>%
  pull(largest_gene) -> largest_genes_to_plot

p3<- qc.metrics %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
  arrange(largest_gene) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) +
  geom_point(size=1) +
  scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1")))
p4<- qc.metrics %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
  arrange(largest_gene) %>%
  ggplot(aes(x=complexity_diff, y=percent.Largest.Gene, colour=largest_gene)) +
  geom_point()+
  scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1")))


write.csv( largest_gene_list, paste0(outdir, sampleid, "largest_geneList.csv"))

png(paste0(outdir, sampleid, ".complexityQC.png"), width = 1152, height = 576)
ggarrange(p1+p2+p3+p4)
dev.off()

###################

#########################################################
### SeuratSct transform pre-proc
#########################################################
print(paste("Seurat initial clustering:",Sys.time()))

adata_filt = adata

adata_filt <- SCTransform(adata_filt, verbose = FALSE) #, vars.to.regress = "percent.mt")
adata_filt <- RunPCA(adata_filt, verbose = FALSE)
adata_filt <- RunUMAP(adata_filt, dims = 1:50, verbose = FALSE)
adata_filt <- FindNeighbors(adata_filt, dims = 1:50, verbose = FALSE)
adata_filt <- FindClusters(adata_filt, algorithm=1, resolution = .3, verbose=FALSE)

p1 <- DimPlot(adata_filt, reduction='umap', group.by='seurat_clusters', label=TRUE, label.size=8, repel=TRUE) 
p1 <- p1 + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('RNA / SCT transf')


#p2 <- VlnPlot(adata_filt, features='nCount_SCT', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata_filt$nCount_SCT), linetype='dashed') + 
#  theme(plot.title = element_text(size = 15), axis.text = element_text(size=15))
#p3 <- VlnPlot(adata_filt, features='nFeature_SCT', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata_filt$nFeature_SCT), linetype='dashed') + 
#  theme(plot.title = element_text(size = 15), axis.text = element_text(size=15))
p4 <- VlnPlot(adata_filt, features='nCount_RNA', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata_filt$nCount_RNA), linetype='dashed') + 
  theme(plot.title = element_text(size = 15), axis.text = element_text(size=15))
p5 <- VlnPlot(adata_filt, features='nFeature_RNA', group.by='seurat_clusters', pt.size=0, log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata_filt$nFeature_RNA), linetype='dashed') + 
  theme(plot.title = element_text(size = 15), axis.text = element_text(size=15))

### cell cycle scoring
adata_filt<-CellCycleScoring(adata_filt, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

p6<- DimPlot(adata_filt, reduction='umap', group.by='Phase', label=TRUE, label.size=8, repel=TRUE) + xlab('UMAP 1') + ylab('UMAP 2') + ggtitle('Cell cycle')


png(paste0(outdir, sampleid, ".Prelim.clustering.png"), width = 700, height = 576)
ggarrange(p1, p6, p4, p5)
dev.off()

qcplots <-FeaturePlot(adata_filt, features = c('nCount_RNA','nFeature_RNA','percent.mt',
                                     'nCount_SCT','nFeature_SCT','percent.Largest.Gene',
                                     'scDblFinder.score',  'S.Score', "G2M.Score"), cols = c("gray80", "red3") )

markerplots<-FeaturePlot(adata_filt, features = c('HIST1H1B', 'MCM6',  "CDC20", 'GNLY',"IFNG",
                                     "NCAM1","IL32", "GZMB", "KIT" ,"NCR1", "CD3G", "NCAM1"))


png(paste0(outdir, sampleid, ".featurePlotQC.png"), width = 900, height = 576)
qcplots
dev.off()

png(paste0(outdir, sampleid, ".featurePlotGenes.png"), width = 1100, height = 576)
markerplots
dev.off()

######################################
### Second doublet calc using clustered data using Doublet Finder
######################################

print(paste("Doublet Finder analysis:",Sys.time()))

## pK Identification (no ground-truth) ------------------------------
sweep.res.list <- paramSweep_v3(adata_filt, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

PK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

homotypic.prop <- modelHomotypic(adata_filt$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(adata_filt@meta.data))                    ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

adata_filt <- doubletFinder_v3(adata_filt, PCs = 1:10, pN = 0.25, pK = PK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)


n = ncol(adata_filt@meta.data)
colnames(adata_filt@meta.data)[c(n-1,n)] = c('pANN_init','DFclass.1')

adata_filt <- doubletFinder_v3(adata_filt, PCs = 1:10, pN = 0.25, pK = 0.08, nExp = nExp_poi.adj, 
                              reuse.pANN = "pANN_init", sct = TRUE)

colnames(adata_filt@meta.data)[n+1] = c('DFclass.2')

adata_filt@meta.data$scDblFinder.class =stringr::str_to_title(adata_filt@meta.data$scDblFinder.class)


png(paste0(outdir, sampleid, ".DoubletFinderDimplot.png"), width = 1000, height = 300)

p <- DimPlot(adata_filt, reduction='umap', 
             group.by=c('DFclass.1', 'DFclass.2', 'scDblFinder.class')
             , label=FALSE, label.size=8, repel=TRUE) 
p <- p + xlab('UMAP 1') + ylab('UMAP 2') 

p
dev.off()

##### COmpare with scdblfinder run before filters
df = data.frame( scDb = adata_filt@meta.data$scDblFinder.class =='Doublet',
                 DF1 = adata_filt@meta.data$DFclass.1 =='Doublet',
                 DF2 = adata_filt@meta.data$DFclass.2 =='Doublet')

vc = limma::vennCounts(df) 

pdf(paste0(outdir, sampleid, ".DoubletFinderOutput.pdf"), width = 8, height = 4)
par(mfrow=c(1,2), mar= c(4,4,2,0))
plot(x = as.numeric(as.character(bcmvn$pK)), main = "Doublet Finder pk",
     y = as.numeric(bcmvn$BCmetric), pch = 16, 
     col = "#41b6c4", xlab="pK", ylab="BCmetric", type="o")

text (y=as.numeric(max(bcmvn$BCmetric)), x= PK, col="red", label=paste("pk=", PK), pos=4)

limma::vennDiagram(vc, cex = 1)
dev.off()

######################################
### SoupX correction
######################################

print(paste("SoupX analysis:",Sys.time()))

DefaultAssay(adata_filt) <- 'RNA'
toc = GetAssayData(object = adata_filt, slot = "counts") ## filtered matrix
tod = GetAssayData(CreateSeuratObject(counts = fdata), slot = "counts") ### unfiltered matrix [ rows were renames ]

sc = SoupChannel(tod,toc)

metadata <- (cbind(as.data.frame(adata_filt[["umap"]]@cell.embeddings),
                   as.data.frame(Idents(adata_filt)),
                   as.data.frame(Idents(adata_filt))))
colnames(metadata) <- c("RD1","RD2","Cluster","Annotation")

sc = setDR(sc,metadata[colnames(sc$toc),c("RD1","RD2")])
sc = setClusters(sc,setNames(metadata$Cluster,rownames(metadata)))

pdf(paste0(outdir, sampleid, ".SoupXRhoEst.pdf"), width = 6, height = 4)

sc = autoEstCont(sc,tfidfMin=0.5)
dev.off()

#check out the genes most expressed in the background
#they suggest not just taking this list, but is interesting to see
write.csv( 
head(sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),],n=20)
          , paste0(outdir, sampleid, "SoupXtop20Bkg_geneList.csv"))

#options(repr.plot.width=10, repr.plot.height=5)
#plotMarkerDistribution(sc)

out = adjustCounts(sc, roundToInt = TRUE)

p1= plotMarkerMap(sc,'GNLY')
p2= plotChangeMap(sc, out, "GNLY")

p3= plotMarkerMap(sc,'GZMB')
p4= plotChangeMap(sc, out, "GZMB")

p5= plotMarkerMap(sc,'MT-CO3')
p6= plotChangeMap(sc, out, "MT-CO3")


png(paste0(outdir, sampleid, ".SoupXExampleCorrection.png"), width = 750, height = 750)
ggarrange(p1+p2, p3+p4,p5+p6, ncol=1)
dev.off()

### Load Post SoupX RNA into Seurat Object 
adata2 = CreateSeuratObject(out,  project =sampleid )
adata2[['percent.mt']] <- PercentageFeatureSet(adata2, pattern = '^MT-')

#add in previous raw RNA data as another assay (RNA_raw)

raw_rna_assay <- CreateAssayObject(counts = toc)
adata2[['RNA_raw']] <- raw_rna_assay


adata2@meta.data <- cbind(adata2@meta.data,adata_filt@meta.data[,c('DFclass.1', 'scDblFinder.class', 'Phase')])

#sum(rownames(adata_filt@meta.data)!=rownames(adata2@meta.data))

write.csv( adata_filt@meta.data , paste0(outdir, sampleid, ".filt.metadata.csv"))

### Save final Seurat object to an .rds file and output the list of keep BCs
print(paste("Saving RDS and barcodes list:",Sys.time()))
saveRDS(adata2, file = paste0(outdir,sampleid, '.intermediate_filtered.rds'))
filtered_bcs <- colnames(adata[["RNA"]])
write(filtered_bcs, file=paste0(outdir,sampleid, '.filtered_barcodes.txt'),sep='\n')
print(paste("Done:",Sys.time()))

write.csv(
data.frame(prefilt= ncol(tod), postFilt = ncol(toc), 
           scDb_doublet= sum(df[,1]), DoubleFinder_doublet=sum(df[,2])),
            paste0(outdir,sampleid, '.stats.csv'))
