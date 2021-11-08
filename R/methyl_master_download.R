#!/usr/bin/env Rscript

####################### For download ####################################
download.barcodes.file="BLCA_tumor_subsample_sample_names.txt"
download.dir="BLCA_idat_files"
download.project=proj
download.data.category="Raw microarray data"
download.data.type="Raw intensities" 
download.exp.strategy="Methylation array" 
download.legacy=TRUE
download.file.type=".idat"
download.platform="Illumina Human Methylation 450"
download.method="api"
download.files.per.chunk=20
download.method="client"
download.barcode <- read.table(paste0(data.dir,
                                      file.sep,
                                      download.barcodes.file))$V1

download.dir <- paste0(data.dir,
                       file.sep,
                       download.folder.name)

download.barcode <- 
  read.table(file = paste0(files.dir,
                           file.sep,
                           download.barcodes.file),
             header=FALSE,
             stringsAsFactors = FALSE)$V1

##also gene and expression data
query <- GDCquery(project       = download.project,
                  data.category          = download.data.category,
                  data.type              = download.data.type,
                  experimental.strategy  = download.exp.strategy, 
                  legacy                 = download.legacy,
                  barcode                = download.barcode,
                  file.type              = download.file.type)

tryCatch(
  GDCdownload(query, 
              method = download.method, 
              files.per.chunk = download.files.per.chunk,
              directory = download.dir),
  error = function(e){
    GDCdownload(query, 
                method = download.method)
  })

if(os.type=="windows"){
  ##For Windows OS:
  shell(paste0('FOR /R \"',
               idat.dir,
               '\" %i IN (*.idat) DO MOVE ',
               '\"%i\" ',
               '\"',
               idat.pooled.files.dir,
               '\"'))
}else if(os.type=="mac" | os.type=="linux"){
  ##For Mac and Linux:
  system(paste0('for i in ',
                idat.dir,
                '/*.idat; ',
                'do mv ',
                '$i ',
                idat.pooled.files.dir))
}else{
  stop("Error: please choose a valid <os.type>")
}
