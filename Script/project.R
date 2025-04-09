#1- Read study design 
#2-creat path abundance file 
#3- download the annotation file 
#4-Tximport to read files

#1
studydesign_read <- read_tsv("studydesign.txt")
sampless <- studydesign_read$sample  
#2
abundance.path <- file.path(sampless) 
#3
GTF <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name")) 
GTF<- as_tibble(GTF)
GTF <- dplyr::rename(GTF, target_id = tx_id) 
GTF<- dplyr::select(GTF, "target_id", "gene_name")
#4
GTF.TXImport <- tximport(abundance.path, 
                     type = "kallisto", 
                     tx2gene = GTF, 
                     txOut = FALSE,
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
myTPM <- GTF.TXImport$abundance
myCounts <- GTF.TXImport$counts
####################################
###data befor filtred and non-normalized

#1
myDGEList <- DGEList(myCounts)
#2
cpm <- cpm(myDGEList)  
#3
log2.cpm <- cpm(myDGEList, log=TRUE) 
#4
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampless) 
#5
tidtable1 <- pivot_longer(log2.cpm.df, #
                                  cols = "WT_FT-194_rep1":"OV151", 
                                  names_to = "samples", 
                                  values_to = "expression")

#6
chart1=ggplot(tidtable1) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95,   
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
####################################
###2- Data filtered and non-normalized

#0
keepers <- rowSums(cpm>1)>=6
myDGEList.keepers <- myDGEList[keepers,]
dim(myDGEList.keepers) 

#1
log2.cpm.2 <- cpm(myDGEList.keepers, log=TRUE)
#2
logafter <- as_tibble(log2.cpm.2, rownames = "geneID")
colnames(logafter) <- c("geneID", sampless)
#
tidtable2 <- pivot_longer(logafter, 
                               cols ="WT_FT-194_rep1":"OV151",  
                                 names_to = "samples", )
                                  values_to = "expression" 
#4
p2=ggplot(tidtable2) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


#################################
#############3- Data filtered and normalized

#1
myDGEList.norm <- calcNormFactors(myDGEList.keepers, method = "TMM")

#2
myDGEList.norm.cpm3 <- cpm(myDGEList.norm, log=TRUE)
#3
log2.cpm.filtered.norm.df <- as_tibble(lmyDGEList.norm.cpm3, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampless)
#4
tidtable3 <- pivot_longer(log2.cpm.filtered.norm.df, 
                                                cols ="WT_FT-194_rep1":"OV151",
                                                names_to = "samples", 
                                                values_to = "expression")

#5
p3=ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

#6
plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)



##############################################################################
###############DEGS volcane polt
#1
group <- factor(studydesign_read$group)
#2
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
#3
cont.matrix <- makeContrasts(health.dis=disease-healthy,levels=design)
#4
v <- voom(myDGEList.filtered.norm,design,plot = TRUE)
#5
fit <- lmFit(v)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
write.fit(fit.cont, file="lmfit_results.txt") 
dim(v$E)
#6
myTopHits <- topTable(fit.cont, adjust ="BH", coef=1, number=11329, sort.by="logFC") 
myTopHits.df <- myTopHits %>% 
  as_tibble(rownames = "geneID")
#7
ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, ) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) + 
  geom_vline(xintercept = 2, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -2, linetype="longdash", colour="#2C467A", size=1) +
  #annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="Red") +
  #annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="Blue") +
  labs(title="Volcano plot",
       subtitle = "Cutaneous leishmaniasis",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

####################################################### 
###############heatmap
#1
col_rs <- bluered(75)
#2
results <- decideTests(fit.cont, method="global", adjust.method="BH", p.value=0.01, lfc=2)
colnames(v$E) <- sampless
diffGenes <- v$E[results[,1] !=0,]
dim(diffGenes)                  
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") 
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") 
module.assign <- cutree(clustRows, k=2)
#3
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=rev(col_rs), scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20)) 

###func
1#
myTopHits <- topTable(fit.cont, adjust ="BH", coef=1, number=50, sort.by="logFC")
gost.res <- gost(rownames(myTopHits), organism = "hsapiens", correction_method = "fdr")
data=gost.res$result
2#
gostplot(gost.res, interactive = T, capped = F)
mygostplot=gostplot(gost.res, interactive = F, capped = F)
3#
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0071944", "GO:0010033"),
  filename = NULL,
  width = NA,
  height = NA)
