library(magrittr)
# group -------------------------------------------------------------------
geo_group <- function(pdata, sample_id, group_id, treat_name, con_name, treat_key, con_key){
  geogroup <- pdata %>% 
    dplyr::select(ID = sample_id, group0 = group_id) %>%
    dplyr::mutate(group = dplyr::case_when(
      grepl(treat_key, group0) ~ treat_name,
      grepl(con_key, group0) ~ con_name
    )) %>% 
    dplyr::select(1,3) %>% 
    dplyr::filter(group == treat_name | group == con_name) %>% 
    dplyr::arrange(group)
  if(geogroup$ group[1] != con_name){
    geogroup <- geogroup %>% 
      dplyr::arrange(desc(group))
  }
  message("Please check the order of samples!")
  return(geogroup)
}

# gpl ---------------------------------------------------------------------
gpl_pro <- function(gpl_name, probe_id, gene_id)
{
  gpl_name <- gpl_name %>% 
    dplyr::select(probe_id, gene_id) %>% 
    dplyr::filter(!grepl("///", .[,2])) %>% 
    dplyr::filter(.[,2] != "") %>% 
    dplyr::filter(!grepl("---", .[,2])) %>% 
    dplyr::mutate(ID = as.character(ID)) %>% 
    dplyr::rename(GENE_SYMBOL = gene_id) %>%  
    dplyr::select(probe_id, GENE_SYMBOL)
  return(gpl_name)
}


# geo_matrix --------------------------------------------------------------
geo_matrix <- function(data_matrix, data_group, gpl_name)
  {
  data_matrix <- data_matrix[, data_group$ID] %>% 
    data.frame() %>% 
    dplyr::mutate(across(!where(is.numeric), as.numeric)) %>% 
    dplyr::mutate(ID = row.names(.) %>% as.character()) %>% 
    dplyr::left_join(., gpl_name, by = "ID") %>% 
    na.omit()
  gene_row <- data_matrix$GENE_SYMBOL
  data_matrix <- data_matrix %>% 
    dplyr::select(-ID, -GENE_SYMBOL) %>% 
    as.matrix()
  rownames(data_matrix) <- gene_row
  data_matrix <- data_matrix %>%  
    limma::avereps() %>%
    data.frame()
  print(dim(data_matrix))
  return(data_matrix)
}

# log2 --------------------------------------------------------------------

geo_log2 <- function(data)
{
  print(c("before:",range(data)))
  ex <- as.data.frame(data)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  
  if (LogC) {
    for(i in 1:ncol(ex)){
      ex[which(ex[,i] < 0),i] <- NaN
    }
    data <- log2(data + 1) 
    print("log2 transform finished")
    print(c("after:",range(data)))
  }else{
    print("log2 transform not needed")
  }
  return(data)
}

# mytheme -----------------------------------------------------------------
library(ggplot2)
mytheme <- 
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        plot.title = element_text(size = 7)) + 
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(),
        panel.borde = element_rect(fill = NA, size = 0.75 * 0.47)) + 
  theme(axis.line = element_line(size = 0.75 * 0.47), 
        axis.text = element_text(size = 6, color = "black"), 
        axis.title = element_text(size = 6),
        axis.ticks = element_line(size = 0.75 * 0.47)) +
  theme(legend.key = element_rect(fill = "white"),
        legend.key.size = unit(c(0.3, 0.3), "cm"), 
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(),
        legend.box.margin = margin(),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_blank(), legend.spacing = unit(0, "cm"),
        legend.box.background = element_blank())
