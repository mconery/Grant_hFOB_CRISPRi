module load R/3.3.2
R

library(tidyverse)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(gplots)
library(ggrepel)

dir="/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/DE/hFOB"
setwd(dir)
htseq_dir = "/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/HTSeq/"

#### READ COUNT ######
geneCount=read.table(file.path(htseq_dir,"htseq_count.txt"), sep="\t", stringsAsFactors=F,row.names=1,header=T)
geneRemove=unlist(tibble(gene=rownames(geneCount)) %>% filter(grepl("__",gene)),use.names=F)
geneCount=geneCount[!rownames(geneCount) %in% geneRemove,]

#get hFOBs
geneCount = geneCount[,grepl("hFOB", names(geneCount))]

#Library size before removing rNRNA
apply(geneCount,2,sum)

# remove rRNA 
gene_id_rm = bind_rows(
        read_delim(
                "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.rRNA_gene_id.txt",
                delim="\t", col_names=c("gene_id","type")
        ),
        read_delim(
                "/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.chrM_gene_id.txt",
                delim="\t", col_names=c("gene_id","type")
        )
) %>% distinct(gene_id)
# remove chrM and rRNA
geneCount = geneCount[!(rownames(geneCount) %in% gene_id_rm$gene_id),]

apply(geneCount,2,sum)


geneCount =  geneCount[,grepl("_N", names(geneCount))==F] #N samples too low, Jim said they were repeated for that reason
names(geneCount) = gsub("hFOB_diff", "hFOBdiff",names(geneCount))
names(geneCount) = gsub("_Q", "", names(geneCount))

###
## After here is only relevant to differential gene expression besides saving geneCount
###

##### DESIGN #######
samples=colnames(geneCount)
exp_design = tibble(sample=samples) %>% separate(sample, c("condition","rep"),remove=F) # change
group=factor(exp_design$condition)
group=relevel(group,ref="hFOBs")
######### FILTER LOW EXPRESSION GENES ############
y=DGEList(counts=geneCount,group=group)
y_cpm = cpm(y)
dim(y$samples) #6 2
dim(y$counts)  #  50054     6
range(y$samples$lib.size)  # 25,936,411 76,417,326
top_sample_perc=(dim(y$samples)[1]-dim(y$samples)[2]+1)/dim(y$samples)[1]
 ### assume top 3 out of 12 are from the same condition
quantile_cpm=apply(y_cpm,1,function(x){quantile(x,top_sample_perc)})
quantile_cpm_quantile=quantile(
        quantile_cpm,
        c(0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8)
)
quantile_cpm_quantile

cpm_cutoff= 1.103826906 #### filter out lowest 50% ~60% genes
range(y$samples$lib.size)/10^6*cpm_cutoff  # assume you need as little as ~ min_count reads in smallest samples to prove that gene express
# 49.05256 66.0390

# filter now
keep=(quantile_cpm > cpm_cutoff)
y <- y[keep, , keep.lib.sizes=FALSE]
range(y$samples$lib.size)  # 27682100 94827229 (27M to 94.8M)
dim(y$counts) # 15547 genes   9 samples

######### calculate normalization factor (TMM) ############
y <- calcNormFactors(y,method="TMM")


################# QC plots ###################
y_cpm=data.frame(gene=rownames(cpm(y)),cpm(y)) %>% as_tibble  #use normalized lib.size read_count/(norm.factor*lib.size)

y_cpm_log=data.frame(gene=rownames(cpm(y,log=T)),cpm(y,log=T)) %>% as_tibble

y_cpmVSN=normalizeVSN(y)
y_cpmVSN=data.frame(gene=rownames(y_cpmVSN),y_cpmVSN) %>% as_tibble

y_expAll_melt=bind_rows(
        y_cpm %>% gather(sample,exp, -gene) %>% mutate(exp_type="cpm"),
        y_cpm_log %>% gather(sample,exp, -gene) %>% mutate(exp_type="cpm_log"),
        y_cpmVSN %>% gather(sample,exp, -gene) %>% mutate(exp_type="cpmVSN")
) %>% arrange(gene)
y_expAll_melt=y_expAll_melt %>% spread(exp_type,exp)

### MSDplot
library(plyr)
sd2mean <- ddply(y_expAll_melt, .(gene), summarize,
cpm_mean=mean(cpm),
cpm_sd=sd(cpm),
lcpm_mean=mean(cpm_log),
lcpm_sd=sd(cpm_log),
cpmVSN_mean=mean(cpmVSN),
cpmVSN_sd=sd(cpmVSN)
)
detach(package:plyr)

sd2mean <- sd2mean %>% tbl_df %>% gather(type,value, -gene) %>% separate(type, c("valueType","xy"), sep="_")
sd2mean <- sd2mean %>% spread(xy, value)

save(geneCount,y,y_expAll_melt,sd2mean,exp_design,file=file.path(dir,"data.RData"))

# plot sdVSmean
p1 <- ggplot(sd2mean,aes(mean,sd)) +
facet_wrap(~valueType, scale="free") +
geom_point() +
geom_smooth(color='red')
theme(axis.text.x = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"))
dir.create(file.path(dir, "qcPlots"), showWarnings = FALSE)
ggsave(file.path(dir,"qcPlots/y_plotSDmean.pdf"), plot=p1, width = 21, height = 7)


################################### PCA plot ################################
myPCAPlot <- function(df, type) {
        myvars <- apply(as.matrix(df %>% select(gene, sample, (!!type)) %>% spread(sample, (!!type)) %>% select(-gene)), 1, var)
        d <- df %>% select(gene, sample, (!!type)) %>% spread(sample, (!!type)) %>% mutate(sd=myvars) %>% filter(sd!=0) %>% select(-gene, -sd)
        pca_m <- prcomp(t(d),
                         center = TRUE,
                         scale. = TRUE)
        pca_df <- data.frame(pca_m$x, reps=row.names(pca_m$x), stringsAsFactors = F) %>% tbl_df
        pca_df <- pca_df %>% separate(reps,c("condition","rep_num"),sep="_", remove=F)
        eigs <- pca_m$sdev^2
        pca_summary <- cbind(
          SD = sqrt(eigs),
          Proportion = eigs/sum(eigs),
          Cumulative = cumsum(eigs)/sum(eigs)) %>% tbl_df
        pca_summary <- pca_summary %>%
          mutate(component=row_number())
        PC1_lab=unlist(pca_summary %>% 
                filter(component %in% c(1,2,3)) %>% 
                mutate(name=paste0("PC",component," (",sprintf("%1.2f%%", 100*Proportion),")")) %>% 
                filter(component==1) %>% 
                select(name),use.names=F
        )
        PC2_lab=unlist(pca_summary %>% 
                        filter(component %in% c(1,2,3)) %>% 
                        mutate(name=paste0("PC",component," (",sprintf("%1.2f%%", 100*Proportion),")")) %>% 
                        filter(component==2) %>% 
                        select(name),use.names=F
                )
        PC3_lab=unlist(pca_summary %>% 
                        filter(component %in% c(1,2,3)) %>% 
                        mutate(name=paste0("PC",component," (",sprintf("%1.2f%%", 100*Proportion),")")) %>% 
                        filter(component==3) %>% 
                        select(name),use.names=F
                )
                p1 <- ggplot(pca_df, aes(PC1,PC2)) +
                  geom_point(aes(col=condition), size=4) +
                  theme_minimal()+
                  theme(axis.text.y   = element_text(size=14),
                        axis.text.x   = element_text(size=14),
                        axis.title.y  = element_text(size=14),
                        axis.title.x  = element_text(size=14),
                        panel.border = element_rect(colour = "black", fill=NA, size=2))+
                  xlab(PC1_lab) + ylab(PC2_lab) +
                  geom_text_repel(aes(label=rep_num))
                #   geom_text(aes(label=reps), check_overlap = TRUE, size = 2)
                p2 <- ggplot(pca_df, aes(PC1,PC3)) +
                  geom_point(aes(col=condition), size=4) +
                  xlab(PC1_lab) + ylab(PC3_lab) +
                  theme_minimal()+
                  theme(axis.text.y   = element_text(size=14),
                        axis.text.x   = element_text(size=14),
                        axis.title.y  = element_text(size=14),
                        axis.title.x  = element_text(size=14),
                        panel.border = element_rect(colour = "black", fill=NA, size=2))+
                  geom_text_repel(aes(label=rep_num))
                p3 <- ggplot(pca_df, aes(PC2,PC3)) +
                    geom_point(aes(col=condition), size=4) +
                    xlab(PC2_lab) + ylab(PC3_lab) +
                    theme_minimal()+
                    theme(axis.text.y   = element_text(size=14),
                          axis.text.x   = element_text(size=14),
                          axis.title.y  = element_text(size=14),
                          axis.title.x  = element_text(size=14),
                          panel.border = element_rect(colour = "black", fill=NA, size=2))+
                    geom_text_repel(aes(label=rep_num))
                p4 <- ggplot(pca_summary) +
                  geom_histogram(aes(component,Cumulative),
                                 stat="identity", fill="yellow", color="yellow", alpha=0.5) +
                  geom_hline(yintercept = 0.8, color="red") +
                  ylab('Cumulative proportion of variance')
                # grid.arrange(p1,p2,ncol=1)
                list(p1,p2,p3,p4,pca_df,pca_summary,PC1_lab,PC2_lab,PC3_lab)
}

##### all samples
p <- myPCAPlot(y_expAll_melt, "cpmVSN")
pdf(file.path(dir,"qcPlots/all_PCA_3.pdf"),width=21, height=7)
do.call("grid.arrange", c(p[1:3], ncol=3))
dev.off()


# ### heatmap
myHeatmap <- function (df) {
  sample_cor <- as.matrix(cor(df))
  heatmap.2(sample_cor,
            distfun = function(x) as.dist((1-x)/2),
            hclust=function(x) hclust(x,method="complete"),
            Rowv = TRUE, Colv= T,
            col=brewer.pal(9,"Blues"),
            symm = T,
            margins = c(12,12),
            trace="none",
            density.info="none",
            key.title = NA,
            keysize=1,
            key.xlab="correlation")
}

mysds <- apply(as.matrix(y_expAll_melt %>% select(gene, sample, cpmVSN) %>% spread(sample, cpmVSN) %>% select(-gene)), 1, sd)
d <- y_expAll_melt %>% select(gene, sample, cpmVSN) %>% spread(sample, cpmVSN) %>% mutate(sd=mysds) %>% filter(sd!=0) %>% select(-gene, -sd)

pdf(file.path(dir,'qcPlots/y_cpmVSN_correlation_heatmap.pdf'))
myHeatmap(d)
dev.off()



d <-  semi_join(
        y_expAll_melt %>% select(gene, sample, cpmVSN) %>% spread(sample, cpmVSN) %>% mutate(sd=mysds) %>% filter(sd!=0),
        sd2mean %>% filter(valueType=="cpmVSN") %>% arrange(desc(mean)) %>% head(n=500)
) %>% select(-gene, -sd)
pdf(file.path(dir,'qcPlots/y_cpmVSN_correlation_heatmap_top500.pdf'))
myHeatmap(d)
dev.off()




#Differential analysis 

## estimateDisp and glmQLFit
### design matrix
groups <- factor(exp_design$condition)
design <- model.matrix(~0+groups)
colnames(design) = gsub("groups", "", colnames(design))
colnames(design)
y <- estimateDisp(y, design) 
fit <- glmQLFit(y, design, robust=T)

my.contrasts <- makeContrasts(
    hFob_restrict.vs.hFob_perm = hFOBdiff-hFOBs,
    levels=design
)

## run glmQLFTest in loop
writeTest <- function(name,dir){
        results = glmQLFTest(fit,contrast=my.contrasts[,name])
        FDR = p.adjust(results$table$PValue, 'fdr')
        df = data.frame(gene_id=rownames(results$table),results$table, FDR=FDR) %>% tbl_df %>% arrange(FDR)
        filename = paste0(name,".DE_edgeR.txt")
        write.table(df, file=file.path(dir,filename), sep="\t", row.name=F, quote=F)
        df %>% mutate(comp=name)
}

dir.create(file.path(dir,"tables"),showWarnings = FALSE)
DE_list = lapply(colnames(my.contrasts),writeTest,dir=file.path(dir,"tables"))
DE_df = do.call("bind_rows",DE_list)
out = DE_df %>% 
filter(FDR < 0.01 & abs(logFC)>0.5849625) %>% 
mutate(direction=ifelse(logFC > 0,"up","down")) %>% 
group_by(direction,comp) %>% summarise(DEgene_num=n_distinct(gene_id)) %>% 
spread(key=direction,value=DEgene_num) %>% mutate(total=down+up)
out
 # comp               down    up total
 # <chr>             <int> <int> <int>
# hFob_restrict.vs.hFob_perm     5413  5734 11147

write.table(out, file=file.path(dir,"tables","pairwise_comparison_DEgeneNum_summary_samplesRemoved.txt"), 
sep="\t", quote=F, row.names=F)


#Get annotation
gene_len = read.delim("/mnt/isilon/sfgi/suc1/customerized_geneome_annotation/hg19/genecode_v19/gencode_name.txt", stringsAsFactors=FALSE, header=F)
names(gene_len) = c("gene_id", "gene_name")
hg19_name_key = gene_len[,c(1,2)]

DElist.n = lapply(DE_list, function(x){
    left_join(DE_df, hg19_name_key)
  })

DElist.n = do.call("rbind",DElist.n)
DElist.n[complete.cases( DElist.n)==F,]$gene_name = DElist.n[complete.cases( DElist.n)==F,]$gene_name = DElist.n[complete.cases( DElist.n)==F,]$gene_name = DElist.n[complete.cases( DElist.n)==F,]$gene_id
write.csv(DElist.n , file="tables/hFOB_DifferentialAnalysis_RNASeq.csv", quote=F, row.names=F)



