# FROM ubuntu:latest

# apt-get update && apt-get install -y

# RUN apt install r-cran-devtools

# FROM rocker/r-ver:4.1.2 

# FROM rocker/r-base:latest

FROM rocker/tidyverse:latest

# FROM rocker/r-ubuntu:20.04

COPY /source_packages /packages

# FROM mcr.microsoft.com/windows/servercore:ltsc2019 as microsoft

# FROM mcr.microsoft.com/windows/nanoserver:1809

# rocker/r-base:latest

ARG R_VERSION

ENV R_VERSION=${R_VERSION}

# RUN curl.exe -o Rtools.exe https://cran.r-project.org/bin/windows/Rtools/rtools40-x86_64.exe

# RUN Rtools.exe /TYPE=full /VERYSILENT -NoNewWindow -Wait

# RUN wget https://centos.pkgs.org/7/springdale-computational-x86_64/R-core-4.1.2-1.sdl7.x86_64.rpm.html

# FROM nuest/rocker-win:ltsc2019-latest

# RUN Rscript -e 'install.packages("devtools", repos = "http://cran.us.r-project.org")'

# WORKDIR .

# RUN install2.r --error devtools

# RUN install2.r --error remotes

# WORKDIR /root/packages/AnnotationDbi

# RUN Rscript -e 'print(getwd()); setwd(getwd()); devtools::install(pkg="/root/packages/AnnotationDbi",            type="source", dependencies=TRUE, upgrade="never")'

# RUN install2.r --error -r http://bioconductor.org/packages/3.14/bioc --deps TRUE AnnotationDbi \
#    && rm -rf /tmp/downloaded_packages/

#RUN install.r "packages/AnnotationDbi/AnnotationDbi_1.56.2.tar.gz"

#RUN build.r "packages/AnnotationDbi/AnnotationDbi_1.56.2.tar.gz"

#RUN Rscript -e 'if(file.exists("packages/htmlwidgets_1.5.4.tar.gz")){utils::install.packages("packages/htmlwidgets_1.5.4.tar.gz", repos = NULL, type="source")}else{stop("problem")}'

#RUN Rscript -e 'print(getwd()); devtools::install("packages/Biobase",                         type="source", dependencies=TRUE, upgrade="never")'


Run Rscript -e 'BiocManager::install("data.table")'
RUN Rscript -e 'devtools::install_local("packages/BiocGenerics_0.36.1.tar.gz")'                            
RUN Rscript -e 'devtools::install_local("packages/BiocParallel_1.24.1.tar.gz")'                                                                                                
RUN Rscript -e 'devtools::install_local("packages/cowplot_1.1.1.tar.gz")'                                
# RUN Rscript -e 'devtools::install_local("packages/data.table_1.14.2.zip")'
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
# apt-get update
# apt-get install -y libbz2-dev
# apt-get install -y liblzma-dev
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update && \
apt-get install -y libbz2-dev && \
apt-get install -y liblzma-dev && \
rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'BiocManager::install("Rhtslib")'
RUN Rscript -e 'BiocManager::install("S4Vectors")'
RUN Rscript -e 'devtools::install_local("packages/Esean-cho-Epicopy-1.08-17-g0eda5f8.tar.gz")'                                            
# RUN Rscript -e 'devtools::install_github("sean-cho/Epicopy")'
# RUN Rscript -e 'devtools::install_local("packages/EpicopyData_1.0.1.tar.gz")'  
Run Rscript -e 'install.packages("packages/EpicopyData_1.0.1.tar.gz", repos = NULL, type = "source")'
# Epicopy installs minfi
# RUN Rscript -e 'devtools::install_local("packages/TCGAbiolinks_2.22.4.tar.gz")'
RUN Rscript -e 'BiocManager::install("TCGAbiolinks")'
# RUN Rscript -e 'devtools::install_local("packages/ungeviz_0.1.0.tar.gz")'
RUN Rscript -e 'devtools::install_github("clauswilke/ungeviz")'
# RUN Rscript -e 'devtools::install_local("packages/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0.tar.gz")'   
Run Rscript -e 'BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")'        
# RUN Rscript -e 'devtools::install_local("packages/CNVRanger_1.6.1.tar.gz")'
RUN Rscript -e 'BiocManager::install("CNVRanger")'
# Above installs RaggedExperiment
# RUN Rscript -e 'devtools::install_local("packages/mknoll-cnAnalysis450k-a60de17.tar.gz")'
RUN Rscript -e 'devtools::install_github("mknoll/cnAnalysis450k")'
# RUN Rscript -e 'devtools::install_local("packages/sdchandra-CNAclinic-3ce2099.tar.gz")'                                    
RUN Rscript -e 'devtools::install_github("sdchandra/CNAclinic")'
# RUN Rscript -e 'devtools::install_local("packages/ENmix_1.26.10.tar.gz")' 
RUN Rscript -e 'BiocManager::install("ENmix")'   
# Run Rscript -e 'devtools::install_local("NourMarzouka-CopyNumber450kCancer-5f24c41.tar.gz")'  
Run Rscript -e 'devtools::install_github("NourMarzouka/CopyNumber450kCancer")'                              
# RUN Rscript -e 'devtools::install_local("packages/matter_1.16.0.tar.gz")'
RUN Rscript -e 'BiocManager::install("matter")'                                                      
# RUN Rscript -e 'devtools::install_local("packages/ChAMP_2.20.0.tar.gz")'
RUN Rscript -e 'BiocManager::install("ChAMP")'
RUN Rscript -e 'devtools::install_local("packages/ChAMPdata_2.26.0.tar.gz")'
RUN Rscript -e 'devtools::install_local("packages/sesame_1.10.5.tar.gz")'
# RUN Rscript -e 'BiocManager::install("sesame", version=3.13)'
RUN Rscript -e 'devtools::install_local("packages/sesameData_1.12.0.tar.gz")'
RUN Rscript -e 'devtools::install_github("sean-cho/Epicopy")'
RUN Rscript -e 'BiocManager::install("peakRAM")'
RUN Rscript -e 'BiocManager::install("ExperimentHub")'
Run Rscript -e 'devtools::install_github("bmbolstad/preprocessCore", dependencies = T, upgrade = "always", configure.args = "--disable-threading")'
RUN Rscript -e 'devtools::install_github("mmariani123/MethylMasteR")'
                   
# When I run and try to install methylmaster there i get: ERROR: dependencies 
# ‘ChAMP’, ‘org.Hs.eg.db’, ‘Epicopy’, ‘TCGAbiolinks’, ‘sesame’, ‘minfi’, ‘CNVRanger’, 
# ‘cnAnalysis450k’, ‘CNAclinic’, ‘IlluminaHumanMethylation450kanno.ilmn12.hg19’, 
# ‘ENmix’, ‘ungeviz’, ‘ChAMPdata’, ‘CopyNumber450kCancer’, ‘RaggedExperiment’, 
# ‘matter’ are not available for package ‘MethylMasteR’

# I tried again and now I have:
# ERROR: dependencies ‘Epicopy’, ‘TCGAbiolinks’, ‘CNVRanger’, ‘cnAnalysis450k’, 
# ‘CNAclinic’, ‘ENmix’, ‘RaggedExperiment’

# docker run --rm -p 127.0.0.1:8787:8787 -e USER=mike -e PASSWORD=mike -e ROOT=TRUE methylmaster

# RUN R -e "install.packages('BiocManager')"
# RUN R -e "install.packages('devtools')"
# RUN R -e "BiocManager::install('AnnotationDbi',                               update=FALSE)"
# RUN R -e "BiocManager::install('Biobase',                                     update=FALSE)"
# RUN R -e "BiocManager::install('BiocGenerics',                                update=FALSE)"
# RUN R -e "BiocManager::install('BiocParallel',                                update=FALSE)"
# RUN R -e "BiocManager::install('ChAMP',                                       update=FALSE)"
# RUN R -e "BiocManager::install('ChAMPdata',                                   update=FALSE)"
# RUN R -e "BiocManager::install('CNAclinic',                                   update=FALSE)"
# RUN R -e "BiocManager::install('cnAnalysis450k',                              update=FALSE)"
# RUN R -e "BiocManager::install('CNVRanger',                                   update=FALSE)"
# RUN R -e "BiocManager::install('cowplot',                                     update=FALSE)"
# RUN R -e "BiocManager::install('data.table',                                  update=FALSE)"
# RUN R -e "BiocManager::install('DelayedArray',                                update=FALSE)"
# RUN R -e "BiocManager::install('dendextend',                                  update=FALSE)"
# RUN R -e "BiocManager::install('DNAcopy',                                     update=FALSE)"
# RUN R -e "BiocManager::install('doParallel',                                  update=FALSE)"
# RUN R -e "BiocManager::install('dplyr',                                       update=FALSE)"
# RUN R -e "BiocManager::install('ENmix',                                       update=FALSE)"
# RUN R -e "BiocManager::install('Epicopy',                                     update=FALSE)"
# RUN R -e "BiocManager::install('EpicopyData',                                 update=FALSE)"
# RUN R -e "BiocManager::install('ExperimentHub',                               update=FALSE)"
# RUN R -e "BiocManager::install('foreach',                                     update=FALSE)"
# RUN R -e "BiocManager::install('future',                                      update=FALSE)"
# RUN R -e "BiocManager::install('GenomicRanges',                               update=FALSE)"
# RUN R -e "BiocManager::install('GEOquery',                                    update=FALSE)"
# RUN R -e "BiocManager::install('ggplot2',                                     update=FALSE)"
# RUN R -e "BiocManager::install('gtools',                                      update=FALSE)"
# RUN R -e "BiocManager::install('htmlwidgets',                                 update=FALSE)"
# RUN R -e "BiocManager::install('igraph',                                      update=FALSE)"
# RUN R -e "BiocManager::install('IlluminaHumanMethylation450kanno.ilmn12.hg19',update=FALSE)"
# RUN R -e "BiocManager::install('illuminaio',                                  update=FALSE)"
# RUN R -e "BiocManager::install('iterators',                                   update=FALSE)"
# RUN R -e "BiocManager::install('limma',                                       update=FALSE)"
# RUN R -e "BiocManager::install('limma',                                       update=FALSE)"
# RUN R -e "BiocManager::install('magrittr',                                    update=FALSE)"
# RUN R -e "BiocManager::install('marray',                                      update=FALSE)"
# RUN R -e "BiocManager::install('Matrix',                                      update=FALSE)"
# RUN R -e "BiocManager::install('matter',                                      update=FALSE)"
# RUN R -e "BiocManager::install('minfi',                                       update=FALSE)"
# RUN R -e "BiocManager::install('org.Hs.eg.db',                                update=FALSE)"
# RUN R -e "BiocManager::install('pheatmap',                                    update=FALSE)"
# RUN R -e "BiocManager::install('plyr',                                        update=FALSE)"
# RUN R -e "BiocManager::install('profmem',                                     update=FALSE)"
# RUN R -e "BiocManager::install('profvis',                                     update=FALSE)"
# RUN R -e "BiocManager::install('RaggedExperiment',                            update=FALSE)"
# RUN R -e "BiocManager::install('ramify',                                      update=FALSE)"
# RUN R -e "BiocManager::install('readxl',                                      update=FALSE)"
# RUN R -e "BiocManager::install('rlist',                                       update=FALSE)"
# RUN R -e "BiocManager::install('S4Vectors',                                   update=FALSE)"
# RUN R -e "BiocManager::install('scales',                                      update=FALSE)"
# RUN R -e "BiocManager::install('sesame',                                      update=FALSE)"
# RUN R -e "BiocManager::install('sesameData',                                  update=FALSE)"
# RUN R -e "BiocManager::install('shape',                                       update=FALSE)"
# RUN R -e "BiocManager::install('TCGAbiolinks',                                update=FALSE)"
# RUN R -e "BiocManager::install('ungeviz',                                     update=FALSE)"
# RUN R -e "devtools::install_github('mmariani123/methylmaster', upgrade="never")"

# WORKDIR ../

# RUN R CMD INSTALL packages/AnnotationDbi_1.56.2.tar.gz                               --type=source
# RUN R CMD INSTALL packages/Biobase_2.50.0.tar.gz                                     --type=source
# RUN R CMD INSTALL packages/BiocGenerics_0.36.1.tar.gz                                --type=source
# RUN R CMD INSTALL packages/BiocParallel_1.24.1.tar.gz                                --type=source
# RUN R CMD INSTALL packages/ChAMP_2.20.0.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/ChAMPdata_2.26.0.tar.gz                                   --type=source
# RUN R CMD INSTALL packages/CNAclinic_1.0.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/cnAnalysis450k_0.99.26.tar.gz                             --type=source
# RUN R CMD INSTALL packages/CNVRanger_1.6.1.tar.gz                                    --type=source
# RUN R CMD INSTALL packages/cowplot_1.1.1.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/data.table_1.14.2.zip                                     --type=source
# RUN R CMD INSTALL packages/DelayedArray_0.16.3.tar.gz                                --type=source
# RUN R CMD INSTALL packages/dendextend_1.15.1.tar.gz                                  --type=source
# RUN R CMD INSTALL packages/DNAcopy_1.68.0.tar.gz                                     --type=source
# RUN R CMD INSTALL packages/doParallel_1.0.16.tar.gz                                  --type=source
# RUN R CMD INSTALL packages/dplyr_1.0.7.tar.gz                                        --type=source
# RUN R CMD INSTALL packages/ENmix_1.26.10.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/Epicopy.tar.gz                                            --type=source
# RUN R CMD INSTALL packages/EpicopyData_1.0.1.tar.gz                                  --type=source
# RUN R CMD INSTALL packages/ExperimentHub_1.16.0.tar.gz                               --type=source
# RUN R CMD INSTALL packages/ExperimentHub_1.16.1.tar.gz                               --type=source
# RUN R CMD INSTALL packages/foreach_1.5.1.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/future_1.21.0.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/GenomicRanges_1.42.0.tar.gz                               --type=source
# RUN R CMD INSTALL packages/GEOquery_2.62.1.tar.gz                                    --type=source
# RUN R CMD INSTALL packages/ggplot2_3.3.3.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/gtools_3.8.2.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/htmlwidgets_1.5.4.tar.gz                                  --type=source
# RUN R CMD INSTALL packages/igraph_1.2.6.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0.tar.gz --type=source
# RUN R CMD INSTALL packages/illuminaio_0.32.0.tar.gz                                  --type=source
# RUN R CMD INSTALL packages/iterators_1.0.13.tar.gz                                   --type=source
# RUN R CMD INSTALL packages/limma_3.30.13.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/limma_3.50.0.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/magrittr_2.0.1.tar.gz                                     --type=source
# RUN R CMD INSTALL packages/marray_1.72.0.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/Matrix_1.3-2.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/matter_1.16.0.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/minfi_1.40.0.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/org.Hs.eg.db_3.14.0.tar.gz                                --type=source
# RUN R CMD INSTALL packages/pheatmap_1.0.12.tar.gz                                    --type=source
# RUN R CMD INSTALL packages/plyr_1.8.6.tar.gz                                         --type=source
# RUN R CMD INSTALL packages/profmem_0.6.0.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/profvis_0.3.7.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/RaggedExperiment_1.14.2.tar.gz                            --type=source
# RUN R CMD INSTALL packages/ramify_0.3.3.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/readxl_1.3.1.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/rlist_0.4.6.2.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/S4Vectors_0.28.1.tar.gz                                   --type=source
# RUN R CMD INSTALL packages/scales_1.1.1.tar.gz                                       --type=source
# RUN R CMD INSTALL packages/sesame_1.8.10.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/sesame_1.10.4.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/sesame_1.11.15.tar.gz                                     --type=source
# RUN R CMD INSTALL packages/sesameData_1.12.0.tar.gz                                  --type=source
# RUN R CMD INSTALL packages/shape_1.4.6.tar.gz                                        --type=source
# RUN R CMD INSTALL packages/TCGAbiolinks_2.22.4.tar.gz                                --type=source
# RUN R CMD INSTALL packages/ungeviz_0.1.0.tar.gz                                      --type=source
# RUN R CMD INSTALL packages/mmariani123-MethylMasteR-ac220fd.tar.gz                   --type=source

# The below opens up R session each time it looks like:

# RUN R.exe --no-save 'install.packages("packages/AnnotationDbi_1.56.2.tar.gz",                               repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/Biobase_2.50.0.tar.gz",                                     repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/BiocGenerics_0.36.1.tar.gz",                                repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/BiocParallel_1.24.1.tar.gz",                                repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ChAMP_2.20.0.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ChAMPdata_2.26.0.tar.gz",                                   repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/CNAclinic_1.0.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/cnAnalysis450k_0.99.26.tar.gz",                             repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/CNVRanger_1.6.1.tar.gz",                                    repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/cowplot_1.1.1.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/data.table_1.14.2.zip",                                     repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/DelayedArray_0.16.3.tar.gz",                                repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/dendextend_1.15.1.tar.gz",                                  repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/DNAcopy_1.68.0.tar.gz",                                     repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/doParallel_1.0.16.tar.gz",                                  repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/dplyr_1.0.7.tar.gz",                                        repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ENmix_1.26.10.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/Epicopy.tar.gz",                                            repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/EpicopyData_1.0.1.tar.gz",                                  repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ExperimentHub_1.16.0.tar.gz",                               repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ExperimentHub_1.16.1.tar.gz",                               repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/foreach_1.5.1.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/future_1.21.0.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/GenomicRanges_1.42.0.tar.gz",                               repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/GEOquery_2.62.1.tar.gz",                                    repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ggplot2_3.3.3.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/gtools_3.8.2.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/htmlwidgets_1.5.4.tar.gz",                                  repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/igraph_1.2.6.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.0.tar.gz", repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/illuminaio_0.32.0.tar.gz",                                  repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/iterators_1.0.13.tar.gz",                                   repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/limma_3.30.13.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/limma_3.50.0.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/magrittr_2.0.1.tar.gz",                                     repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/marray_1.72.0.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/Matrix_1.3-2.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/matter_1.16.0.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/minfi_1.40.0.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/org.Hs.eg.db_3.14.0.tar.gz",                                repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/pheatmap_1.0.12.tar.gz",                                    repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/plyr_1.8.6.tar.gz",                                         repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/profmem_0.6.0.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/profvis_0.3.7.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/RaggedExperiment_1.14.2.tar.gz",                            repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ramify_0.3.3.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/readxl_1.3.1.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/rlist_0.4.6.2.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/S4Vectors_0.28.1.tar.gz",                                   repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/scales_1.1.1.tar.gz",                                       repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/sesame_1.8.10.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/sesame_1.10.4.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/sesame_1.11.15.tar.gz",                                     repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/sesameData_1.12.0.tar.gz",                                  repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/shape_1.4.6.tar.gz",                                        repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/TCGAbiolinks_2.22.4.tar.gz",                                repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/ungeviz_0.1.0.tar.gz",                                      repos = NULL, type="source")'
# RUN R.exe --no-save 'install.packages("packages/mmariani123-MethylMasteR-ac220fd.tar.gz",                   repos = NULL, type="source")'
