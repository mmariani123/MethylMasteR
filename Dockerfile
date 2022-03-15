FROM rocker/tidyverse:latest

COPY /source_packages /packages

ARG R_VERSION

ENV R_VERSION=${R_VERSION}

Run Rscript -e 'BiocManager::install("data.table")'
RUN Rscript -e 'devtools::install_local("packages/BiocGenerics_0.36.1.tar.gz")'                            
RUN Rscript -e 'devtools::install_local("packages/BiocParallel_1.24.1.tar.gz")'                                                                                                
RUN Rscript -e 'devtools::install_local("packages/cowplot_1.1.1.tar.gz")'                                
RUN Rscript -e 'devtools::install_local("packages/plyr_1.8.6.tar.gz")'                               
RUN Rscript -e 'devtools::install_local("packages/dplyr_1.0.7.tar.gz")'
RUN Rscript -e 'devtools::install_local("packages/magrittr_2.0.1.tar.gz")'                                   
RUN Rscript -e 'devtools::install_local("packages/DelayedArray_0.16.3.tar.gz")'                                
RUN Rscript -e 'devtools::install_local("packages/dendextend_1.15.1.tar.gz")'                                                                   
RUN Rscript -e 'devtools::install_local("packages/doParallel_1.0.16.tar.gz")'                                
RUN Rscript -e 'devtools::install_local("packages/ExperimentHub_1.16.0.tar.gz")'                                                           
RUN Rscript -e 'devtools::install_local("packages/foreach_1.5.1.tar.gz")'                                      
RUN Rscript -e 'devtools::install_local("packages/future_1.21.0.tar.gz")'                                      
RUN Rscript -e 'devtools::install_local("packages/GenomicRanges_1.42.0.tar.gz")'                               
RUN Rscript -e 'devtools::install_local("packages/GEOquery_2.62.1.tar.gz")'                                    
RUN Rscript -e 'devtools::install_local("packages/ggplot2_3.3.3.tar.gz")'                                      
RUN Rscript -e 'devtools::install_local("packages/gtools_3.8.2.tar.gz")'                                       
RUN Rscript -e 'devtools::install_local("packages/htmlwidgets_1.5.4.tar.gz")'                                  
RUN Rscript -e 'devtools::install_local("packages/igraph_1.2.6.tar.gz")'                                       
RUN Rscript -e 'devtools::install_local("packages/illuminaio_0.32.0.tar.gz")'                                  
RUN Rscript -e 'devtools::install_local("packages/iterators_1.0.13.tar.gz")'                                                                      
RUN Rscript -e 'devtools::install_local("packages/limma_3.50.0.tar.gz")'                                                                        
RUN Rscript -e 'devtools::install_local("packages/marray_1.72.0.tar.gz")'                                      
RUN Rscript -e 'devtools::install_local("packages/Matrix_1.3-2.tar.gz")'                                                                                                                                           
RUN Rscript -e 'devtools::install_local("packages/pheatmap_1.0.12.tar.gz")'                                                                         
RUN Rscript -e 'devtools::install_local("packages/profmem_0.6.0.tar.gz")'                                      
RUN Rscript -e 'devtools::install_local("packages/profvis_0.3.7.tar.gz")'                                                                
RUN Rscript -e 'devtools::install_local("packages/ramify_0.3.3.tar.gz")'                                       
RUN Rscript -e 'devtools::install_local("packages/readxl_1.3.1.tar.gz")'                                       
RUN Rscript -e 'devtools::install_local("packages/rlist_0.4.6.2.tar.gz")'                                                                       
RUN Rscript -e 'devtools::install_local("packages/scales_1.1.1.tar.gz")'
RUN Rscript -e 'devtools::install_local("packages/shape_1.4.6.tar.gz")'                                                                    
RUN Rscript -e 'devtools::install_local("packages/DNAcopy_1.68.0.tar.gz")'                                                                           
RUN Rscript -e 'devtools::install_local("packages/org.Hs.eg.db_3.14.0.tar.gz")' 
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && \
apt-get install -y libbz2-dev && \
apt-get install -y liblzma-dev && \
rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'BiocManager::install("Rhtslib")'
RUN Rscript -e 'BiocManager::install("S4Vectors")'
RUN Rscript -e 'devtools::install_local("packages/Esean-cho-Epicopy-1.08-17-g0eda5f8.tar.gz")'                                             
Run Rscript -e 'install.packages("packages/EpicopyData_1.0.1.tar.gz", repos = NULL, type = "source")'
RUN Rscript -e 'BiocManager::install("TCGAbiolinks")'
RUN Rscript -e 'devtools::install_github("clauswilke/ungeviz")' 
Run Rscript -e 'BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")'        
RUN Rscript -e 'BiocManager::install("CNVRanger")'
RUN Rscript -e 'devtools::install_github("mknoll/cnAnalysis450k")'                                   
RUN Rscript -e 'devtools::install_github("sdchandra/CNAclinic")'
RUN Rscript -e 'BiocManager::install("ENmix")'   
Run Rscript -e 'devtools::install_github("NourMarzouka/CopyNumber450kCancer")'                              
RUN Rscript -e 'BiocManager::install("matter")'                                                      
RUN Rscript -e 'BiocManager::install("ChAMP")'
RUN Rscript -e 'devtools::install_local("packages/ChAMPdata_2.26.0.tar.gz")'
RUN Rscript -e 'devtools::install_local("packages/sesame_1.10.5.tar.gz")'
RUN Rscript -e 'devtools::install_local("packages/sesameData_1.12.0.tar.gz")'
RUN Rscript -e 'devtools::install_github("sean-cho/Epicopy")'
RUN Rscript -e 'BiocManager::install("peakRAM")'
RUN Rscript -e 'BiocManager::install("ExperimentHub")'
Run Rscript -e 'devtools::install_github("bmbolstad/preprocessCore", dependencies = T, upgrade = "always", configure.args = "--disable-threading")'
RUN Rscript -e 'devtools::install_github("mmariani123/MethylMasteR")'
