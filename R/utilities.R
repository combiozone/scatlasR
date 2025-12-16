

#' CompareMarkers
#'
#' @param data object
#' @param col_ct cell column name
#' @param col_gene gene column name
#'
#' @return data frame
#' @export
#'
CompareMarkersByPublicDB <- function(
    data=NULL, 
    fdb=NULL,
    col_ct='cellName',
    col_gene='geneSymbol',
    n = 10
) { 
    d_db <- readxl::read_excel(fdb)
    
    df <- NULL
    for (x in unique(data$cluster)){
        print(x)
        geneset <- data$gene[data$cluster == x]
        print(length(geneset))

        for (ct in unique(d_db[[col_ct]])){
            markers <- d_db[[col_gene]][ d_db[[col_ct]] == ct ]
            markers <- str_split(markers, ',')[[1]]
            overlap_genes <- intersect(geneset, markers)
            geneSymbol <- paste(overlap_genes, collapse = ",")
            avg_log2fc <- round(mean(data$avg_log2FC[data$cluster == x & data$gene %in% overlap_genes]), 3)
            n <- length(overlap_genes)
            if (n == 0){ next }
            new_row <- c(cluster=x, cellType=ct, geneCount=n, avg_log2fc=avg_log2fc, geneSymbol=geneSymbol)

            df <- rbind(df, new_row)
        }
    }

    df <- as.data.frame(df)
    print(dim(df))

    # show top10
    d_topX <- df %>% group_by(cluster) %>% arrange(desc(as.numeric(geneCount))) %>% slice_head(n = n)

    return(d_topX)
}







