library(tidyverse)

#' Title
#'
#' @param barcode_dir Where the barcode files for each well are as .txt
#' Some might be missing due to nothing found earlier in the pipeline. 
#' @param output_name Name of the file. It will be created in the current directory.
#' @return
#' @export
#'
#' @examples
aggregate_barcodes <- function(barcode_dir="sample_barcodes_files_rev_primer",
                               output_name="rev_primer_bc_matrix.tsv"){
  
  files <- list.files(barcode_dir, full.names = T) %>% grep("\\.txt",., value = T)
  
  # plate_well names
  # pattern based on GP Edyta's sequencing
  pos_names <- files %>% gsub("^.*(S000.*)_HL.*","\\1",.) %>%
    gsub("[[:alnum:]]+_BC[[:alnum:]]+_","",.)
  
  all_files_list <- lapply(files, readLines)
  # names not kept
  names(all_files_list) <-pos_names
  
  bc_uniques <- all_files_list %>% flatten %>% unique
  
  bc_mat <- matrix(0, nrow=length(bc_uniques), ncol=length(all_files_list))
  rownames(bc_mat) <- bc_uniques
  colnames(bc_mat) <- names(all_files_list)
  
  for(cell in names(all_files_list)){
    ctab <- all_files_list[[cell]] %>% table
    bc_mat[names(ctab),cell] <- bc_mat[names(ctab),cell]+ctab
  }
  
  # edyta barcodes are 52 or 53 length (53 got trimmed earlier, I think)
  bc_mat_only_complete <- bc_mat[rownames(bc_mat) %>% nchar %>% `>=`(52), ]
  
  write.table(bc_mat, output_name, sep="\t")
  write.table(bc_mat_only_complete, gsub("(\\.[[:alpha:]]+)","_only_complete\\1",output_name), sep="\t") 
}

if(sys.nframe()==0){
  args=commandArgs(trailingOnly=T)
  if(length(args)==2){
    aggregate_barcodes(args[1],args[2])
  }
  else {
    message(">> Aggregating barcodes error: Missing arguments! <<")
    
  }
}
