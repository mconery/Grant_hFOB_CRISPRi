library(tidyverse)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(ggrepel)
library(venn)
library(tidyverse)


dir="/mnt/isilon/sfgi/pahlm/analyses/grant/rnaSeq/bone_diff/DE/hFOB"
setwd(dir)


load(file.path(dir,"data.RData")) # save(geneCount,y,y_rmNaive,y_expAll_melt,sd2mean,exp_design,exp_design_rmNaive,file="data.RData")

geneCount_df=data.frame(gene_id=rownames(geneCount), geneCount, stringsAsFactors=F) %>% tbl_df

geneCount_df = geneCount_df %>% gather(key="sample",value="rawCount",-gene_id)

#Use Chun's effective gene length file
gene_len=read_delim("/mnt/isilon/sfgi/programs/HTSeq-0.6.1/geneModels/gencodeV19.lincRNA.snomiRNA.annotation_for_HTseq.length.txt",delim="\t",col_names=c("gene_id","gene_name","gene_len"))
geneCount_df = left_join(geneCount_df,gene_len)

###### calculate TPM
# step 1: Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK)
geneCount_df = geneCount_df %>% mutate(RPK=rawCount/gene_len*1000)

# step 2: Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor
scale_factor=geneCount_df %>% group_by(sample) %>% summarise(total_RPK=sum(RPK)) %>% ungroup() %>% mutate(scale_factor=total_RPK/1000000)

# step3: Divide the RPK values by the “per million” scaling factor. This gives you TPM.
tpm=left_join(geneCount_df, scale_factor) %>% mutate(tpm=RPK/scale_factor) %>% select(gene_id,gene_name,sample,tpm)

# step4 (optional): calculate log2tpm using pseudo-count 1 added to raw count
log2tpm=left_join(geneCount_df, scale_factor) %>% 
mutate(RPK_p=(rawCount+1)/gene_len*1000) %>% 
mutate(log2tpm=log2(RPK_p/scale_factor)) %>% 
select(gene_id,gene_name,sample,log2tpm)

tpm=left_join(tpm,log2tpm)


save(tpm, geneID2name, file=file.path(dir,"tpm_data.Rdata"))
