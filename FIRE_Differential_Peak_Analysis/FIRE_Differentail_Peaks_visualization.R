setwd(getwd()) 
source("./RNA_Seq_Ara/theme_linhua.R")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hiAnnotator))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(hiAnnotator))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL) {
  library(matrixStats)
  library(rtracklayer)
  
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  
  ################
  # Get the column means
  ################
  cm = rowMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = rowSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}

##-------------------------------------------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-f", "--file"), type = "character", help = "Path to first BED file (A)"),
  make_option(c("-c", "--cluster"), type = "integer", default = 6,
              help = "Number of clusters expected [default: %default]"),
  make_option(c("-s", "--order"), type = "character", default = "1,2,3,4,5,6",
              help = "Cluster sort, e.g. \"1,2,3,4,5,6\" [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = "Heatmap.pdf",
              help = "Output PDF filename [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))


if (is.null(opt$file)) {
  stop("Required arguments missing! Usage:
       Rscript script.R -f <file.bed> -c <number> -s \"3,4,1,2,5\" [-o output.pdf]", 
       call. = FALSE)
}


##-------------------------------------------------------------------------------------------------------------------

file_name <- opt$file
ID <- file_name %>%basename() %>% str_remove("\\.[^.]*$")

prefix <- sub("\\.(bed|txt)$", "", file_name, ignore.case = TRUE)

col_count <- ncol(fread(file_name, nrows = 0))

Fiber_raw <- fread(file_name, select = c(1, 2:col_count), header = TRUE) 

Fiber_raw$ID <- paste0("peak", 1:nrow(Fiber_raw))
Fiber_raw <- as.data.frame(Fiber_raw)
rownames(Fiber_raw) <- Fiber_raw$ID
Fiber_raw$ID <- NULL

colnames(Fiber_raw) <- colnames(Fiber_raw) %>%
  basename() %>%
  str_replace_all(pattern = ".bw", replacement = "") %>%
  str_replace_all(pattern = "Sam_", replacement = "")%>%
  str_replace_all(pattern = "'", replacement = "")%>%
  str_replace_all(pattern = "#chr", replacement = "seqnames")

count_data <- round(as.matrix(Fiber_raw %>%as.data.frame() %>% select(4:col_count)))
count_data[is.na(count_data)] <- 0

ExCPM <- edgeR::cpm(count_data)

colnames(ExCPM) <- paste0(colnames(ExCPM), ".cpm")



prefixes <- c() 
ExCPM <- as.data.frame(ExCPM)
for (i in seq(1, col_count-3, by = 2)) {
  if (i + 1 <= col_count) {
    col1 <- colnames(ExCPM)[i]
    col2 <- colnames(ExCPM)[i + 1]
    
    max_common <- ""
    for (j in seq_len(min(nchar(col1), nchar(col2)))) {
      if (substr(col1, j, j) == substr(col2, j, j)) {
        max_common <- substr(col1, 1, j)
      } else {
        break
      }
    }
    

    if (grepl("-", max_common)) {
      dash_pos <- max(gregexpr("-", max_common)[[1]])
      common_prefix <- substr(max_common, 1, dash_pos)
    } else {
      common_prefix <- max_common
    }
    
    prefixes <- c(prefixes, common_prefix)
    prefixes <- prefixes[grepl("-", prefixes)]
    avg_col <- rowMeans(ExCPM[, c(i, i+1)])
    ExCPM[[paste0(common_prefix, "mean.cpm")]] <- avg_col
  }
}



mean_col <- colnames(ExCPM %>%select(contains("mean")))

Fiber_log2 <- ExCPM %>%
  mutate(across(
    c(mean_col),
    ~ log2(. + 1),
    .names = "{.col}_log2"
  ))


mean_log2_col <- colnames(Fiber_log2 %>% select(contains("log2")))

#input <- paste(mean_log2_col, collapse = ", ")
#quoted_input <- paste0('"', strsplit(input, ",\\s*")[[1]], '"', collapse = ", ")
#temp <- rowScale(as.matrix(Fiber_log2[ , ..quoted_input, with = FALSE]))
Fiber_zscore <- as.matrix(Fiber_log2[, mean_log2_col])%>% rowScale()
colnames(Fiber_zscore) <- gsub("\\.cpm_log2$", ".zscore", colnames(Fiber_zscore))

sign <- gsub("-mean.zscore", "", colnames(Fiber_zscore))

setnames(as.data.frame(Fiber_zscore), colnames(Fiber_zscore), sign)

#Fiber_zscore <- as.data.frame(Fiber_zscore)
#Fiber_zscore$ID <- paste0("peak", 1:nrow(Fiber_zscore))
#Fiber_zscore <- as.data.frame(Fiber_zscore)
#rownames(Fiber_zscore) <- Fiber_zscore$ID
#Fiber_zscore$ID <- NULL

mat <- as.matrix(Fiber_zscore)

set.seed(123)
kclust <- kmeans(mat, centers = opt$cluster)
cluster_assign <- data.frame(ID = rownames(mat), Cluster = kclust$cluster)

pdf_file_name <- paste0(prefix, "_Peak_Z-score-Heatmap-hclust.pdf")

pdf(file = pdf_file_name)

order_levels <- strsplit(opt$order, ",")[[1]]
kclust$cluster <- factor(kclust$cluster, levels = order_levels)

Heatmap(
  mat,
  name = "Z-score",
  cluster_rows = FALSE,
  row_split = kclust$cluster,
  col = colorRamp2(c(-1.5,-0.75, 0, 0.75, 1.5), c("blue","skyblue", "white", "pink","red")),
  show_row_names = FALSE,
  #  right_annotation = ha,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  width = unit(6, "cm")
)

dev.off()



boxpolt_df <- as.data.frame(mat) %>%
  tibble::rownames_to_column(var = "ID") %>%  
  left_join(cluster_assign, by = "ID") 


for (cluster_num in 1:opt$cluster) {
  cluster_rows <- boxpolt_df[boxpolt_df$Cluster == cluster_num, ]
  rownames(cluster_rows) <- cluster_rows$ID
  cluster_rows$ID <- NULL
  
  MCHH <- reshape2::melt(
    as.matrix(cluster_rows[, -ncol(cluster_rows)]),
    varnames = c("ID", "variable"),
    value.name = "value"
  ) %>% 
    as.data.table()
  
  MCHH_WML_nONA <- MCHH[!is.na(value), ][, value := as.numeric(value)]
  
  
  median_data <- MCHH_WML_nONA %>% 
    group_by(variable) %>% 
    summarise(median_value = median(value, na.rm = TRUE))
  
  
  
  pdf(file = paste0(ID, "_Cluster_", cluster_num, "_Peak_Z-score_boxplot.pdf"))
  
  p <- ggboxplot(
    MCHH_WML_nONA,
    x = "variable",
    y = "value",
    color = "variable",
    fill = "variable", 
    alpha = 0.8,      
    palette = c("#00AFBB", "#E7B800", "#FC4E07", "#2E8B57", "#4682B4"),  
    xlab = "Treatment Group",
    ylab = "Z-score",
    legend = "none",
    outlier.shape = NA,  
    width = 0.7    
  ) +
    stat_summary(
      fun = median,
      geom = "point",
      shape = 21,
      size = 5,       
      color = "black", 
      fill = "white",    
      stroke = 1     
    ) +
    geom_line(
      data = median_data,
      aes(x = variable, y = median_value, group = 1), 
      color = "black",
      linewidth = 0.5,
      linetype = "solid",
      alpha = 0.8
    ) +
    stat_boxplot(
      geom = "errorbar",
      width = 0.25,      
      linewidth = 0.5,  
      color = "black"   
    ) +
  scale_y_continuous(
    breaks = seq(-1.5, 1.5, 0.5),
    limits = c(-1.5, 1.5),
    expand = c(0, 0)
  ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8,face = "italic"),
      axis.text.y = element_text(size = 8,face = "italic"),
      axis.title.y = element_text(size = 8,face = "italic"),
      axis.title.x = element_blank(),
      axis.ticks.length = unit(0.25, "cm"),
      legend.text = element_text(size = 8,face = "italic"),
      axis.text = element_text(size = 8,face = "italic"),
      aspect.ratio = 0.618/1, 
      axis.line.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.ticks.x = element_line(color = "black", linewidth = 1),
      axis.ticks.length.y = unit(0.2, "cm"), 
      axis.ticks.length.x = unit(0.2, "cm"), 
      axis.ticks.y = element_line(color = "black", linewidth = 1)
    ) +
    
    labs(title = paste("Cluster", cluster_num, "Expression Z-scores"))
  
  print(p) 
  dev.off()
 
  Cluster_region <- MCHH_WML_nONA %>%
    left_join(
      Fiber_raw %>% 
        tibble::rownames_to_column("ID") %>%  
        select(ID, seqnames, start, end),
      by = "ID"  
    )
  
  write.table(
    Cluster_region,
    file = paste0(ID, "_Cluster_", cluster_num, ".bed"),
    sep = "\t",         
    row.names = FALSE,   
    quote = FALSE     
  )
}

##---------------------------------------------------------------------------------------------------------------------------

