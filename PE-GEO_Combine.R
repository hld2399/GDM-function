rm(list = ls()) 

# options()$repos 
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options()$BioC_mirror
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
# options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

## preparation----

dir.create("raw") 
dir.create("output") 
source("function.R") 


library(GEOquery)
gset1 <- getGEO('GSE75010', destdir = "raw", 
                AnnotGPL = F, 
                getGPL = F) 
gset2 <- getGEO('GSE60438', destdir = "raw", 
                AnnotGPL = F, 
                getGPL = F) 


## GSET1

library(magrittr)
library(limma)
head(exprs(gset1[[1]])[,1:6]) 

### group

pdata1 <- pData(gset1[[1]]) 
colnames(pdata1) ## geo_accession source_name_ch1
table(pdata1$"characteristics_ch1") ## group: Glomerulus of Metastasis kidney, Glomerulus of control kidney
head(pdata1$geo_accession)
group1 <- geo_group(pdata = pdata1, ## clinical
                    sample_id = "geo_accession", group_id = "characteristics_ch1", 
                        treat_name = "PE", con_name = "Control",  
                    treat_key = "diagnosis: PE", con_key = "non-PE") 
dir.create("output/GSE75010") 
data.table::fwrite(group1, "output/GSE75010/GSE75010_group.csv", row.names = F) 
table(group1$group) 


unique(pdata1$"platform_id") 
# GPL1 <- getGEO(unique(pdata1$"platform_id"), destdir = "raw", AnnotGPL = F) 
# gpl1 <- Table(GPL1) 
gpl1 <- read.table("raw/GPL6244-17930.txt",sep = "\t",header = T,check.names = F)
colnames(gpl1) 


gpl1$GENE_SYMBOL <- strsplit(gpl1$gene_assignment, split = " /// ") %>% lapply(function(x){
  x1 <- x
  if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
  if(length(x1) == 0) x1 <- x[1]
  vec_gene = strsplit(x, split = " // ") %>% lapply(function(x) if(length(x) > 1) return(x[2]) else return(x)) %>%
    unlist() %>% unique()
  vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
  vec_gene <- paste0(vec_gene, collapse = " /// ")
  return(vec_gene)
}) %>% unlist()

gpl1 <- gpl_pro(gpl_name = gpl1, probe_id = "ID", gene_id = "GENE_SYMBOL")



data1 <- exprs(gset1[[1]]) 
data1 <- geo_matrix(data_matrix = data1, data_group = group1, gpl_name = gpl1)
data.table::fwrite(data1, "output/GSE75010/GSE75010_Matrix.csv", row.names = T)



data1 <- geo_log2(data1)



boxplot(data1, outline = F, las = 2) 
data1_norm <- normalizeBetweenArrays(data1) %>% 
  data.frame()
data.table::fwrite(data1_norm, "output/GSE75010/GSE75010_Matrix_norm.csv", row.names = T)
boxplot(data1_norm, outline = F, las = 2) 



head(exprs(gset2[[1]])[,1:6])

### group

pdata2 <- pData(gset2[[1]])
colnames(pdata2)
head(pdata2$geo_accession)
table(pdata2$"source_name_ch1")

group2 <- geo_group(pdata = pdata2, 
                    sample_id = "geo_accession", group_id = "source_name_ch1", 
                    treat_name = "PE", con_name = "Control", 
                    treat_key = "pre-eclamptic_decidua", con_key = "normotensive control_decidua") 
dir.create("output/GSE60438") 
data.table::fwrite(group2,"output/GSE60438/GSE60438_group.csv", row.names = F)
table(group2$group)


unique(pdata2$"platform_id")
# GPL2 <- getGEO(unique(pdata2$"platform_id"), destdir = "raw", AnnotGPL = F) 
# gpl2 <- Table(GPL2) 
gpl2 <- data.table::fread("raw/GPL10558-50081.txt",sep = "\t",header = T,na.strings = "",fill = T)
colnames(gpl2) 

# gpl2$GENE_SYMBOL <- strsplit(gpl2$gene_assignment, split = " /// ") %>% lapply(function(x){
#   x1 <- x
#   if(length(x1) > 1) x1 <- x1[!grepl("^OTTHU", x1)]
#   if(length(x1) == 0) x1 <- x[1]
#   vec_gene = strsplit(x, split = " // ") %>% lapply(function(x) if(length(x) > 1) return(x[2]) else return(x)) %>%
#     unlist() %>% unique()
#   vec_gene <- vec_gene[!grepl("^OTTHU", vec_gene)]
#   vec_gene <- paste0(vec_gene, collapse = " /// ")
#   return(vec_gene)
# }) %>% unlist()

gpl2 <- gpl_pro(gpl_name = gpl2, probe_id = "ID", gene_id = "Symbol")



data2 <- exprs(gset2[[1]])
data2 <- geo_matrix(data_matrix = data2, data_group = group2, gpl_name = gpl2)
data.table::fwrite(data2,"output/GSE60438/GSE60438_Matrix.csv", row.names = T)



data2 <- geo_log2(data2)



boxplot(data2, outline = F, las = 2)
data2_norm <- normalizeBetweenArrays(data2) %>% 
  data.frame()
data.table::fwrite(data2_norm,"output/GSE60438/GSE60438_Matrix_norm.csv", row.names = T)
boxplot(data2_norm, outline = F, las = 2) 




same_gene <- intersect(rownames(data1), rownames(data2))



data1 <- data1[same_gene,]
data2 <- data2[same_gene,]




count_matrix <- cbind(data1, data2) %>% 
  na.omit()
boxplot(count_matrix, outline = F, las = 2) 

batch <- c(rep("1", length(group1$group)),
           rep("2", length(group2$group)))


library(sva)
adjusted_counts <- ComBat(count_matrix, batch = batch) 



adjusted_counts <- normalizeBetweenArrays(adjusted_counts)
boxplot(adjusted_counts, outline = F, las = 2)



adjusted_group <- rbind(group1, group2)
adjusted_group_sort <- adjusted_group[order(adjusted_group$group, decreasing = F),] 
table(adjusted_group_sort$group)



dir.create("output/Merged_Dataset")
write.csv(adjusted_group_sort,"output/Merged_Dataset/PE_Datasets_Group.csv", row.names = F)



adjusted_counts_sort <- adjusted_counts[,adjusted_group_sort$ID]
write.csv(adjusted_counts_sort, "output/Merged_Dataset/PE_Datasets_Matrix.csv")



adjusted_group_Disease <- adjusted_group_sort[which(adjusted_group_sort$group!="Control"),]
adjusted_counts_Disease <- adjusted_counts_sort[, adjusted_group_Disease$ID]
write.csv(adjusted_counts_Disease,"output/Merged_Dataset/PE_Datasets_Disease_Matrix.csv")



A <- "GSE75010"
B <- "GSE60438"


gs <- factor( c(rep(A, length(group1$group)),
                rep(B, length(group2$group))), levels = c(A, B))

group_col <- c("#8deaf0","#ffb3bf",  "#EED2EE") 

library(ggplot2)


dat_before_long <- count_matrix %>% 
  t() %>% 
  data.frame() %>% 
  tibble::rownames_to_column("sample") %>% 
  dplyr::mutate(group = gs) %>% 
  dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>% 
  tidyr::gather(key = geneid, value, - c(sample, group))

p_before <- 
  ggplot(dat_before_long, aes(sample, value)) + 
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + 
  scale_fill_manual(values = group_col,name="PE_Datasets") +
  labs(title = "Before") + 
  mytheme + 
  theme(legend.position = 'top', 
        legend.direction = "horizontal", 
        # legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size = 6),
    # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 6) 
    legend.text = element_text(size = 6,vjust = 0.5),
    legend.title = element_text(size = 6,vjust = 0.9),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank() 
    ) + 
  coord_cartesian(clip = "off")

p_before 

ggsave("output/boxplot_before.pdf", plot = p_before,units = "cm",width = 11,height = 5) 



dat_after_long <- adjusted_counts %>% 
  t() %>% 
  data.frame() %>%
  tibble::rownames_to_column("sample") %>% 
  dplyr::mutate(group = gs) %>% 
  dplyr::mutate(sample = factor(.$sample, levels = .$sample)) %>% 
  tidyr::gather(key = geneid, value, - c(sample, group))

p_after <- 
  ggplot(dat_after_long, aes(sample, value)) + 
  geom_boxplot(aes(fill = group), lwd = 0.1 *0.47, outlier.shape = NA) + 
  scale_fill_manual(values = group_col,name="PE_Datasets") +
  labs(title = "After") + 
  mytheme +
  theme(legend.position = 'top', 
        legend.direction = "horizontal", 
        # legend.title = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  theme(axis.text.y = element_text(size = 6),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 6) 
        legend.text = element_text(size = 6,vjust = 0.5),
        legend.title = element_text(size = 6,vjust = 0.9),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank() 
  ) +   
  coord_cartesian(clip = "off") 

p_after 

ggsave("output/boxplot_after.pdf", plot = p_after, units = "cm",width = 11,height = 5) 



count_matrix_pca <- prcomp(t(count_matrix), scale. = T)
count_matrix_pcs <- data.frame(count_matrix_pca$x, group = gs)

## pca1,pca2
percentage <- round(count_matrix_pca$sdev / sum(count_matrix_pca$sdev) * 100, 2)
percentage <- paste(colnames(count_matrix_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

pca_before <- ggplot(count_matrix_pcs, aes(x = PC1,y = PC2, color = group)) +
  geom_point(aes(shape = group), size = 0.75) +
  stat_ellipse(level = 0.95, show.legend = F, geom = "polygon", aes(fill = group),  alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  scale_color_manual(values = group_col,name="PE_Datasets") +
  scale_fill_manual(values = group_col) +
  xlab(percentage[1]) +ylab(percentage[2]) +
  labs(title = "Before") +
  mytheme +
  theme(panel.grid = element_line(colour = "grey90", size = 0.75 * 0.47, linetype = 1),
        plot.title = element_text(size = 7, hjust = 0, vjust = 0.5, margin = unit(c(0.1,0.1,0.1,0.1), "cm")),
        legend.position = 'top', 
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_text(size = 6,vjust = 0.9), legend.key.size = unit(c(0.15, 0.15), "cm")) + 
  coord_cartesian(clip = "off")+guides(shape="none")

pca_before

ggsave("output/pca_before.pdf", plot = pca_before, units = "cm",width = 6, height = 5)



adjust_matrix_pca <- prcomp(t(adjusted_counts), scale. = T) 
adjust_matrix_pcs <- data.frame(adjust_matrix_pca$x, group = gs)
percentage <- round(adjust_matrix_pca$sdev / sum(adjust_matrix_pca$sdev) * 100, 2)
percentage <- paste(colnames(adjust_matrix_pcs),"(", paste(as.character(percentage), "%", ")", sep=""))

pca_after <- ggplot(adjust_matrix_pcs, aes(x = PC1,y = PC2, color = group)) +
  geom_point(aes(shape = group), size = 0.75) +
  stat_ellipse(level = 0.95, show.legend = F, geom = "polygon", aes(fill = group),  alpha = 0.2) +
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 0.75 * 0.47) +
  scale_color_manual(values = group_col,name="PE_Datasets") +
  scale_fill_manual(values = group_col) +
  xlab(percentage[1]) +ylab(percentage[2]) +
  labs(title = "After") +
  mytheme +
  theme(panel.grid = element_line(colour = "grey90", size = 0.75 * 0.47, linetype = 1),
        plot.title = element_text(size = 7, hjust = 0, vjust = 0.5, margin = unit(c(0.1,0.1,0.1,0.1), "cm")),
        legend.position = 'top', 
        legend.direction = "horizontal",
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.title = element_text(size = 6,vjust = 0.9), legend.key.size = unit(c(0.15, 0.15), "cm")) + 
  coord_cartesian(clip = "off")+guides(shape="none")

pca_after

ggsave("output/pca_after.pdf", plot = pca_after, units = "cm",width = 6, height = 5)


# 
# library(patchwork)
# p <- p_before + p_after + (pca_before | pca_after) + plot_layout(ncol = 1) +
#   plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 12))
# 
# p
# ggsave("output/Figure2.pdf", plot = p, units = "cm",width = 17, height = 21)

## save
# save.image() 
save(data1, data2,  file = ".Rdata")

  
  
