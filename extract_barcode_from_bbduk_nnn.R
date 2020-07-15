library(tidyverse)

extract_from_seqs <- function(prep_seqs_file, output_file,
                              mask_length=21,
                              barcode_length=52,
                              forward_bc="GTAA",
                              rev_bc="CATT"){
  seqs <- readLines(prep_seqs_file)
  mask <- rep("N",mask_length) %>% paste(collapse = "")
  barcodes <- character(0)
  barcode_min <- min(nchar(forward_bc), nchar(rev_bc), floor(barcode_length*0.7))
  barcode_min <- 5
  for(s in seqs){
    if(grepl(mask, s)){
      
      ante_post_mask <- strsplit(s, mask)
      ante <- ante_post_mask[[1]][1]
      post <- ante_post_mask[[1]][2]
      
      # compare to fw
      if(!is.na(post) && nchar(post)>barcode_min){
        post_match <- substr(post, 1, nchar(forward_bc)) %>% 
          strsplit("") %>% .[[1]] %>%
          `==`(strsplit(forward_bc,"")[[1]])
        if(sum(post_match)>=(nchar(forward_bc)-1)){
          barcodes <- c(barcodes, substr(post,1,barcode_length))
          next
        }
      }
      
      # compare to rev
      if(!is.na(ante) && nchar(ante)>barcode_min){
        ante <- ante %>% strsplit("") %>% .[[1]] %>% rev %>% paste(collapse = "")
        rev_match <- substr(ante, 1, nchar(rev_bc)) %>% 
          strsplit("") %>% .[[1]] %>%
          `==`(strsplit(rev_bc,"")[[1]])
        if(sum(rev_match)>=(nchar(forward_bc)-1)){
          barcodes <- c(barcodes, substr(ante,1,barcode_length))
          next
        }
      }
      barcodes <- c(barcodes, "XXX") # mask found, but barcode too short ot not found

    }
  }
  
  writeLines(barcodes, output_file)
  
}

if(sys.nframe()==0){
  args=commandArgs(trailingOnly=T)
  if(length(args)<3){
    extract_from_seqs(args[1],args[2])
  }
  else if(length(args)<4){
    extract_from_seqs(args[1],args[2],
                      mask_length = args[3]
                      )
  
  } else{
    extract_from_seqs(args[1], args[2],
                      args[3], args[4],
                      args[5], args[6])
  }
}