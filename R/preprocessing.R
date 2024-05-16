#' read input data file and store in required format
#' 
#' @export
#' @param mutation_file mutation data file; see inst/extdata/examples/*_snv.csv for examples
#' @param outputDir output directory for saving output data
importFiles <- function(mutation_file, 
                        outputDir=NULL,
                        alt_reads_thresh = 0, # to be tested
                        vaf_thresh = 0 # to be tested
                        ) {
  
  if (is.null(outputDir)) {
    outputDir = getwd()
  }
  
  mutation_data = importMutationFileOnly(mutation_file, alt_reads_thresh, vaf_thresh)
  mutation_data$is_cn <- c(rep(0, nrow(mutation_data$y)))
  
  return(mutation_data)
}

#' import mutation file
importMutationFileOnly <- function(mutation_file, alt_reads_thresh = 0, vaf_thresh = 0) {
  data <- read_csv(mutation_file, show_col_types = FALSE)
  data <- data %>% filter(alt_reads/total_reads>=vaf_thresh, alt_reads>=alt_reads_thresh)
  output_data <- list()
  
  output_data$y <- as.matrix(data[c("mutation", "sample", "alt_reads")] %>% pivot_wider(names_from = sample, values_from = alt_reads, values_fill = 0))
  rownames(output_data$y) <- output_data$y[,'mutation']
  output_data$y <- output_data$y[,-1, drop=FALSE]
  rowname = rownames(output_data$y)
  colname = colnames(output_data$y)
  output_data$y <- matrix(as.numeric(output_data$y), ncol = ncol(output_data$y))
  rownames(output_data$y) = rowname
  colnames(output_data$y) = colname
  
  output_data$MutID <- rowname
  
  output_data$n <- as.matrix(data[c("mutation", "sample", "total_reads")] %>% pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0))
  rownames(output_data$n) <- output_data$n[,'mutation']
  output_data$n <- output_data$n[,-1, drop=FALSE]
  rowname = rownames(output_data$n)
  colname = colnames(output_data$n)
  output_data$n <- matrix(as.numeric(output_data$n), ncol = ncol(output_data$n))
  rownames(output_data$n) = rowname
  colnames(output_data$n) = colname
  
  if (any((output_data$y - output_data$n) > 0)) {
    stop("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
  }
  
  if (any(output_data$n==0)) {
    print("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- round(mean(output_data$n))
  }
  
  output_data$icn <- as.matrix(data[c("mutation", "sample", "tumor_integer_copy_number")] %>% pivot_wider(names_from = sample, values_from = tumor_integer_copy_number, values_fill = 2))
  rownames(output_data$icn) <- output_data$icn[,'mutation']
  output_data$icn <- output_data$icn[,2:ncol(output_data$icn)]
  output_data$icn <- matrix(as.numeric(output_data$icn), ncol=ncol(output_data$y))
  output_data$icn <- apply(output_data$icn, 1, function(x) {
    if (all(x == 2)) {
      return(2)
    } else {
      return(mean(x[x != 2]))
    }
  })

  output_data$cncf <- as.matrix(data[c("mutation", "sample", "cncf")] %>% pivot_wider(names_from = sample, values_from = cncf, values_fill = 0))
  rownames(output_data$cncf) <- output_data$cncf[,'mutation']
  output_data$cncf <- output_data$cncf[,-1,drop=FALSE]
  rowname = rownames(output_data$cncf)
  colname = colnames(output_data$cncf)
  output_data$cncf <- matrix(as.numeric(output_data$cncf), ncol = ncol(output_data$cncf))
  rownames(output_data$cncf) = rowname
  colnames(output_data$cncf) = colname
  
  output_data$S = ncol(output_data$y)
  output_data$I = nrow(output_data$y)
  
  output_data$tcn = output_data$icn * output_data$cncf + 2 * ( 1 - output_data$cncf)
  
  if ("major_integer_copy_number" %in% colnames(data)) {
    output_data$mtp <- as.matrix(data[c("mutation", "sample", "major_integer_copy_number")] %>% pivot_wider(names_from = sample, values_from = major_integer_copy_number, values_fill = 1))
    rownames(output_data$mtp) <- output_data$mtp[,'mutation']
    output_data$mtp <- output_data$mtp[,2:ncol(output_data$mtp)]
    output_data$mtp <- matrix(as.numeric(output_data$mtp), ncol=ncol(output_data$y))
    output_data$mtp <- apply(output_data$mtp, 1, function(x) {
      if (all(x == 1)) {
        return(1)
      } else {
        return(mean(x[x != 1]))
      }
    })
  } else {
    output_data$mtp <- estimateMultiplicityMatrix(output_data)[,1]
  }
  
  if ("purity" %in% colnames(data)) {
    output_data$purity <- as.matrix(data[c("mutation", "sample", "purity")] %>% pivot_wider(names_from = sample, values_from = purity, values_fill = 0))
    rownames(output_data$purity) <- output_data$purity[,'mutation']
    output_data$purity <- output_data$purity[,-1, drop=FALSE]
    rowname = rownames(output_data$purity)
    colname = colnames(output_data$purity)
    output_data$purity <- matrix(as.numeric(output_data$purity), ncol = ncol(output_data$purity))
    rownames(output_data$purity) = rowname
    colnames(output_data$purity) = colname
    output_data$purity = colSums(output_data$purity) / colSums(!!output_data$purity)
  } else {
    output_data$purity = rep(0.8, ncol(output_data$y))
  }
  
  return(output_data)
}
