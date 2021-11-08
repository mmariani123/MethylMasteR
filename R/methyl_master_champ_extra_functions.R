#!/usr/bin/env Rscript

# Extra ChAMP functions

#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @param champ.beta
#' @param champ.pheno
#' @param champ.mds.plot
#' @param champ.density.plot
#' @param champ.dendrogram
#' @param champ.pdf.plot
#' @param champ.rplot
#' @param champ.feature.sel
#' @param champ.results.dir
#' @return A matrix of the infile
#' @export
champ.QC <- function(beta        = NULL,
                     pheno       = NULL,
                     mdsPlot     = NULL,
                     densityPlot = NULL,
                     dendrogram  = NULL,
                     PDFplot     = NULL,
                     Rplot       = NULL,
                     Feature.sel = NULL,
                     resultsDir  = NULL
                     ){
    par(mar=c(1,1,1,1))
    message("[===========================]")
    message("[<<<<< ChAMP.QC START >>>>>>]")
    message("-----------------------------")

    if(!file.exists(resultsDir)){dir.create(resultsDir)}

    message("champ.QC Results will be saved in ", resultsDir)

    message("[QC plots will be proceed with ", dim(beta)[1],
            " probes and ", dim(beta)[2], " samples.]\n")

    if(min(beta, na.rm = TRUE) == 0){
        beta[beta == 0] <- 1e-06
        message("[", length(which(beta == 0)),
                paste0(" Zeros dectect in your ",
                "dataset, will be replaced with 0.000001]\n"))
    }

    if(ncol(beta) != length(pheno)){
        stop(paste0("Dimension of DataSet Samples, ",
                    "pheno and name must be the same. ",
                    "Please check your input."))
        message("<< Prepare Data Over. >>")
    }

    if(mdsPlot){
        if(Rplot){
            mdsPlot(beta,
                    numPositions = 1000,
                    sampGroups = pheno,
                    colnames(beta)
            )
            if(PDFplot){
                pdf(paste(resultsDir,
                          "raw_mdsPlot.pdf",
                          sep = "/"),
                    width = 6,
                    height = 4)
                mdsPlot(beta,
                        numPositions = 1000,
                        sampGroups = pheno,
                        colnames(beta)
                )
                dev.off()
            }
            message("<< plot mdsPlot Done. >>\n")
        }
    }
    if(densityPlot){
        if(Rplot)
            densityPlot(beta,
                        sampGroups = pheno,
                        main = paste("Density plot of raw data (",
                                     nrow(beta),
                                     " probes)",
                                     sep = ""),
                        xlab = "Beta")
        if(PDFplot){
            pdf(paste(resultsDir,
                      "raw_densityPlot.pdf",
                      sep = "/"),
                width = 6,
                height = 4)
            densityPlot(beta,
                        sampGroups = pheno,
                        main = paste("Density plot of raw data (",
                                     nrow(beta), " probes)", sep = ""),
                        xlab = "Beta")
            dev.off()
        }
        message("<< Plot densityPlot Done. >>\n")
    }
    if(dendrogram){
        if(Feature.sel == "None"){
            message(paste0("< Dendrogram Plot Feature Selection Method >: ",
                           "No Selection, directly use all CpGs to calculate",
                           " distance matrix."))
            hc <- hclust(dist(t(beta)))
        }
        else if (Feature.sel == "SVD") {
            message(paste0("< Dendrogram Feature Selection Method >: ",
                           "Use top SVD CpGs to calculate distance matrix."))
            SVD <- svd(beta)
            rmt.o <- EstDimRMT(beta - rowMeans(beta))
            M <- SVD$v[, 1:rmt.o$dim]
            rownames(M) <- colnames(beta)
            colnames(M) <- paste("Component", c(1:rmt.o$dim))
            hc <- hclust(dist(M))
        }
        dend <- as.dendrogram(hc)
        MyColor <- rainbow(length(table(pheno)))
        names(MyColor) <- names(table(pheno))
        labels_colors(dend) <- MyColor[pheno[order.dendrogram(dend)]]
        dend <- dend %>% set("labels_cex", 0.8)
        ##MM change below
        dend <- dend %>%
            set("leaves_pch", 19) %>%
            set("leaves_cex", 0.6) %>%
            set("leaves_col",
                MyColor[pheno[order.dendrogram(dend)]][1])
        if(Rplot){
            plot(dend,
                 center = TRUE,
                 main = paste("All samples before normalization (",
                              nrow(beta), " probes)", sep = ""))
            legend("topright", fill = MyColor, legend = names(MyColor))
        }
        if(PDFplot){
            pdf(paste(resultsDir,
                      "raw_SampleCluster.pdf",
                      sep = "/"),
                width = floor(log(ncol(beta)) * 3),
                height = floor(log(ncol(beta)) * 2)
            )
            plot(dend,
                 center = TRUE,
                 main = paste("All samples before normalization (",
                              nrow(beta), " probes)", sep = "")
            )
            legend("topright", fill = MyColor, legend = names(MyColor))
            dev.off()
        }
        message("<< Plot dendrogram Done. >>\n")
    }
    message("[<<<<<< ChAMP.QC END >>>>>>>]")
    message("[===========================]")
    message("[You may want to process champ.norm() next.]\n")
}

#' GenStatM.R
#'
#' https://rdrr.io/bioc/FEM/src/R/GenStatM.R
#'
#'
#' @param dnaM.m
#' @param pheno.v
#' @param chiptype
#' @return A matrix of the infile
#' @export
GenStatM <- function(dnaM.m,
                     pheno.v,
                     chiptype="450k"
                     ){

    if (chiptype == "450k"){
        data("probe450kfemanno")
        probefemanno <- probe450kfemanno
    }
    else if (chiptype == "EPIC" ){
        data("probeEPICfemanno")
        probefemanno <- probeEPICfemanno
    }
    else{
        print("ERROR: Please indicate correct data type!")
        break
    }

    extractFn <- function(tmp.v, ext.idx) {
        return(tmp.v[ext.idx])
    }
    map.idx <- match(rownames(dnaM.m), probefemanno$probeID);
    probeInfo.lv <- lapply(probefemanno, extractFn, map.idx)
    beta.lm <- list()
    for (g in 1:6) {
        group.idx <- which(probeInfo.lv[[3]] == g)
        tmp.m <- dnaM.m[group.idx, ]
        rownames(tmp.m) <- probeInfo.lv$eid[group.idx];
        sel.idx <- which(is.na(rownames(tmp.m)) == FALSE);
        tmp.m <- tmp.m[sel.idx,];
        nL <- length(factor(rownames(tmp.m)));
        nspg.v <- summary(factor(rownames(tmp.m)),maxsum=nL);
        beta.lm[[g]] <- rowsum(tmp.m,group=rownames(tmp.m))/nspg.v;
        print(paste("Done for regional gene group ", g, sep = ""))
    }
    unqEID.v <- unique(c(rownames(beta.lm[[2]]), rownames(beta.lm[[4]]),
                         rownames(beta.lm[[1]])))
    avbeta.m <- matrix(nrow = length(unqEID.v), ncol = ncol(dnaM.m))
    colnames(avbeta.m) <- colnames(dnaM.m)
    rownames(avbeta.m) <- unqEID.v
    for (gr in c(1, 4, 2)) {
        avbeta.m[match(rownames(beta.lm[[gr]]), rownames(avbeta.m)),
        ] <- beta.lm[[gr]]
    }
    data.m <- avbeta.m
    sampletype.f <- as.factor(pheno.v)
    design.sample <- model.matrix(~0 + sampletype.f)
    colnames(design.sample) <- levels(sampletype.f)
    sampletypes.v <- levels(sampletype.f)
    lmf.o <- lmFit(data.m, design.sample)
    ntypes <- length(levels(sampletype.f))
    ncomp <- 0.5 * ntypes * (ntypes - 1)
    cont.m <- matrix(0, nrow = ncol(design.sample), ncol = ncomp)
    tmp.v <- vector()
    c <- 1
    for (i1 in 1:(ntypes - 1)) {
        for (i2 in (i1 + 1):ntypes) {
            cont.m[i1, c] <- -1
            cont.m[i2, c] <- 1
            tmp.v[c] <- paste(sampletypes.v[i2], "--", sampletypes.v[i1],
                              sep = "")
            c <- c + 1
        }
    }
    rownames(cont.m) <- sampletypes.v
    colnames(cont.m) <- tmp.v
    lmf2.o <- contrasts.fit(lmf.o, cont.m)
    bay.o <- eBayes(lmf2.o)
    top.lm <- list()
    for (c in 1:ncol(cont.m)) {
        top.lm[[c]] <- topTable(bay.o, coef = c, adjust.method = "fdr",
                                number = nrow(data.m))
    }

    return(list(top = top.lm, cont = cont.m, avbeta = avbeta.m))

}

#' DoIntEpi450k
#'
#' ##https://rdrr.io/bioc/FEM/src/R/DoIntEpi450k.R
#'
#'
#' @param statM.o
#' @param adj.m
#' @param c
#' @param dmaMode
#' @return A matrix of the infile
#' @export
DoIntEpi450k <-
    function(statM.o,
             adj.m,
             c,
             dmaMode="avbeta"
             ){

        if(length(grep("[a-zA-Z]",rownames(adj.m)))!=0){
            print("ERROR: The rownames of adj.m should be EntrezID");
            break
        }

        if (dmaMode == "avbeta"){
            avbeta.m <- statM.o$avbeta;
            commonEID.v <- intersect(rownames(adj.m),rownames(avbeta.m));
            mapA.idx <- match(commonEID.v,rownames(adj.m));
            tmpA.m <- adj.m[mapA.idx,mapA.idx];

            mapM.idx <- match(commonEID.v,rownames(avbeta.m));
            tmpM.m <- avbeta.m[mapM.idx,];

            gr.o <- graph.adjacency(tmpA.m,mode="undirected");
            comp.l <- clusters(gr.o);
            ngpc.v <- summary(factor(comp.l$member));
            maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
            maxc.idx <- which(comp.l$member==maxCLid);

            ##get the max connected network
            tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
            gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
            tmpE.m <- get.edgelist(gr.o);
            tmpM.m <- tmpM.m[maxc.idx,];

            #### now extract statistics
            selcol.idx <- match(c("t","P.Value"),colnames(statM.o$top[[c]]));
            ordrow.idx <- match(rownames(tmpM.m),rownames(statM.o$top[[c]]));
            if(length(intersect(c("ID"),colnames(statM.o$top[[c]])))==1){
                ordrow.idx <- match(rownames(tmpM.m),statM.o$top[[c]][,1]);
            }
            statM.m <- statM.o$top[[c]][ordrow.idx,selcol.idx];
            rownames(statM.m) <- rownames(tmpM.m);

            return(list(statM=statM.m,adj=tmpA.m,avbeta=avbeta.m));

        }
        else if(dmaMode == "singleProbe"){
            probeID.v <- statM.o$probeID[[c]];
            commonEID.v <- intersect(rownames(adj.m), names(probeID.v));
            mapA.idx <- match(commonEID.v,rownames(adj.m));
            tmpA.m <- adj.m[mapA.idx,mapA.idx];

            mapM.idx <- match(commonEID.v,names(probeID.v));
            tmpM.v <- probeID.v[mapM.idx];

            gr.o <- graph.adjacency(tmpA.m,mode="undirected");
            comp.l <- clusters(gr.o);
            ngpc.v <- summary(factor(comp.l$member));
            maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
            maxc.idx <- which(comp.l$member==maxCLid);

            #get the max connected network
            tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
            gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
            tmpE.m <- get.edgelist(gr.o);
            tmpM.v <- tmpM.v[maxc.idx];

            #### now extract statistics
            selcol.idx <- match(c("t","P.Value"),colnames(statM.o$top[[c]]));
            ordrow.idx <- match(names(tmpM.v),rownames(statM.o$top[[c]]));
            if(length(intersect(c("ID"),colnames(statM.o$top[[c]])))==1){
                ordrow.idx <- match(names(tmpM.v),statM.o$top[[c]][,1]);
            }
            statM.m <- statM.o$top[[c]][ordrow.idx,selcol.idx];
            rownames(statM.m) <- names(tmpM.v);

            return(list(statM=statM.m,adj=tmpA.m,probeID=probeID.v));

        }
        else{print("Please indicate correct mode for statM.o object!"); break}
    }


#' DoEpiMod
#'
#' https://rdrr.io/bioc/FEM/src/R/DoEpiMod.R
#'
#'
#' @param intEpi.o
#' @param nseeds
#' @param gamma
#' @param nMC=1000
#' @param sizeR.v
#' @param minsizeOUT
#' @param writeOUT
#' @param nameSTUDY
#' @param ew.v
#' @return A matrix of the infile
#' @export
DoEpiMod <-
    function(intEpi.o,
             nseeds=100,
             gamma=0.5,
             nMC=1000,
             sizeR.v=c(1,100),
             minsizeOUT=10,
             writeOUT=TRUE,
             nameSTUDY="X",
             ew.v=NULL){

        PasteVector <- function(v){
            vt <- v[1];
            if(length(v) > 1){
                for(g in 2:length(v)){
                    vt <- paste(vt,v[g],sep=" ")

                }
            }
            vt <- paste(vt," EnD",sep="");
            out.v <- sub(" EnD","",vt);
            out.v <- sub("NA , ","",out.v);
            out.v <- sub(" , NA","",out.v);
            out.v <- sub(" , NA , "," , ",out.v);
            return(out.v);
        }

        Heaviside <- function(v){
            out.v <- v;
            out.v[which(v>=0)] <- 1;
            out.v[which(v<0)] <- 0;
            return(out.v);
        }

        WriteOutPval <- function(pv.v,round.min=3,round.max=50){
            round.v <- round.min:round.max
            th.v <- 10^(-round.v);
            outpv.v <- vector(length=length(pv.v));
            done.idx <- which(pv.v >= th.v[1]);
            outpv.v[done.idx] <- round(pv.v[done.idx],round.min);
            todo.idx <- setdiff(1:length(pv.v),done.idx);
            for(i in todo.idx){
                if(length(which(th.v <= pv.v[i]))>0){
                    outpv.v[i] <-
                      round(pv.v[i],round.v[min(which(th.v <= pv.v[i]))]);
                }
                else{
                    outpv.v[i] <- 0;
                }
            }
            return(outpv.v);
        }

        ##require(igraph);
        ##require(org.Hs.eg.db);

        statM.m <- intEpi.o$statM;
        adj.m <- intEpi.o$adj;
        statM.v <- statM.m[,1];
        nameSTUDY <- paste("Epi-",nameSTUDY,sep="");
        statI.v <- abs(statM.v);
        names(statI.v) <- rownames(statM.m);

        x <- org.Hs.egSYMBOL;
        mapped_genes <- mappedkeys(x)
        xx <- as.list(x[mapped_genes])
        mapEIDtoSYM.v <- unlist(xx);
        map.idx <- match(rownames(adj.m),names(mapEIDtoSYM.v));
        anno.m <- cbind(rownames(adj.m),mapEIDtoSYM.v[map.idx]);
        colnames(anno.m) <- c("EntrezID","Symbol");

        ### find subnetworks around seeds
        ntop <- nseeds;

        ### use greedy approach: rank nodes to define seeds
        rankN.s <- sort(statI.v,decreasing=TRUE,index.return=TRUE);
        seedsN.v <- names(statI.v)[rankN.s$ix[1:ntop]];

        ### now define weights on network
        print("Constructing weighted network");
        tmpA.m <- adj.m;
        gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
        tmpE.m <- get.edgelist(gr.o);
        if(is.null(ew.v)){
            tmpW.v <- vector(length=nrow(tmpE.m));
            for(e in 1:nrow(tmpE.m)){
                match(tmpE.m[e,],rownames(tmpA.m)) -> map.idx;
                tmpW.v[e] <- mean(statI.v[map.idx]);
                print(e);
            }
        }
        else{
            tmpW.v <- ew.v
        }

        ### a number of edges might have a weight of zero which would
        ### later alter the topology of network. this is not desired,
        ### hence we replace 0 by the minimum non-zero value.
        minval <- min(setdiff(tmpW.v,0))
        tmpW.v[tmpW.v==0] <- minval;

        ### defined weighted graph and adjacency matrix
        grW.o <- set.edge.attribute(gr.o,"weight",value=tmpW.v);
        adjW.m <- get.adjacency(grW.o,attr="weight")

        ### Run Spin-Glass algorithm
        print("Running Spin-Glass algorithm");
        sizeN.v <- vector();
        sgcN.lo <- list();
        for(v in 1:ntop){
            sgcN.o <- spinglass.community(gr.o,
                                          weights=tmpW.v,
                                          spins=25,
                                          start.temp=1,
                                          stop.temp=0.1,
                                          cool.fact=0.99,
                                          update.rule=c("config"),
                                          gamma=gamma,
                                          vertex=rankN.s$ix[v]);
            sizeN.v[v] <- length(sgcN.o$comm);
            sgcN.lo[[v]] <- sgcN.o;
            print(paste("Done for seed ",v,sep=""));
        }
        names(sizeN.v) <- seedsN.v;
        print("Module sizes=");
        print(sizeN.v);
        ### compute modularities
        modN.v <- vector();
        for(v in 1:ntop){
            subgr.o <- induced.subgraph(grW.o,sgcN.lo[[v]]$comm);
            modN.v[v] <- mean(get.edge.attribute(subgr.o,name="weight"))
        }
        names(modN.v) <- seedsN.v;
        print("Modularity values=");
        print(modN.v);

        ### now determine significance against randomisation of profiles
        print("Starting Monte Carlo Runs");
        modNmc.m <- matrix(nrow=ntop,ncol=nMC);
        for(m in 1:ntop){
            subgr.o <- induced.subgraph(gr.o,sgcN.lo[[m]]$comm);
            nN <- sizeN.v[m];
            if( (nN> sizeR.v[1]) && (nN< sizeR.v[2])){
                tmpEL.m <- get.edgelist(subgr.o);
                for(run in 1:nMC){
                    permN.idx <-
                      sample(1:nrow(tmpA.m),nrow(tmpA.m),replace=FALSE);
                    tmpEW.v <- vector();
                    for(e in 1:nrow(tmpEL.m)){
                        match(tmpEL.m[e,],rownames(tmpA.m)[permN.idx]) ->
                          map.idx;
                        tmpEW.v[e] <- mean(statI.v[map.idx]);
                    }
                    subgrW.o <-
                      set.edge.attribute(subgr.o,"weight",value=tmpEW.v)
                    modNmc.m[m,run] <-
                      mean(get.edge.attribute(subgrW.o,name="weight"));
                }
            }
            print(paste("Done for seed/module ",m,sep=""));
        }

        modNpv.v <- rep(1,ntop);
        for(v in 1:ntop){
            if( (sizeN.v[v] > sizeR.v[1]) && (sizeN.v[v]< sizeR.v[2])){
                modNpv.v[v] <- length(which(modNmc.m[v,] > modN.v[v]))/nMC;
            }
        }
        names(modNpv.v) <- seedsN.v;
        print(modNpv.v);

        ### summarize hits
        print("Summarising and generating output");
        selpvN.idx <- which(modNpv.v < 0.05);
        selSize.idx <- which(sizeN.v >= minsizeOUT);
        selMod.idx <- intersect(selpvN.idx,selSize.idx);
        print(selMod.idx);
        print(seedsN.v);
        topmodN.m <- matrix(nrow=length(selMod.idx),ncol=6);
        match(seedsN.v[selMod.idx],anno.m[,1]) -> map.idx;
        seedsSYM.v <- anno.m[map.idx,2];

        topmodN.m[,1] <- seedsN.v[selMod.idx];
        topmodN.m[,2] <- seedsSYM.v;
        topmodN.m[,3:5] <- cbind(sizeN.v[selMod.idx],
                                 modN.v[selMod.idx],
                                 modNpv.v[selMod.idx]);
        mi <- 1;
        for(m in selMod.idx){
            tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[m]]$comm];
            genes.v <- anno.m[match(tmpEID.v,anno.m[,1]),2];
            topmodN.m[mi,6] <- PasteVector(genes.v);
            mi <- mi+1;
        }
        colnames(topmodN.m) <- c("EntrezID(Seed)",
                                 "Symbol(Seed)",
                                 "Size",
                                 "Mod",
                                 "P",
                                 "Genes");

        if(writeOUT){
            write.table(topmodN.m,file=paste("topEPI-",
                                             nameSTUDY,
                                             ".txt",
                                             sep=""),
                        quote=FALSE,
                        sep="\t",
                        row.names=FALSE);
        }

        seltopmodN.lm <- list();
        for(m in 1:length(selMod.idx)){
            tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[selMod.idx[m]]]$comm]

            match(tmpEID.v,anno.m[,1]) -> map.idx;

            match(tmpEID.v,rownames(tmpA.m)) -> map1.idx;

            seltopmodN.m <- cbind(anno.m[map.idx,1:2],statM.m[map1.idx,],
                                  statI.v[map1.idx]);

            seltopmodN.lm[[m]] <- seltopmodN.m;

            colnames(seltopmodN.lm[[m]]) <- c("EntrezID",
                                              "Symbol",
                                              "stat(DNAm)",
                                              "P(DNAm)",
                                              "stat(Int)");
        }

        names(seltopmodN.lm) <- seedsSYM.v

        if(writeOUT){

            for(m in 1:length(selMod.idx)){

                out.m <- seltopmodN.lm[[m]];

                out.m[,3] <- round(as.numeric(out.m[,3]),2);

                out.m[,4] <- WriteOutPval(as.numeric(out.m[,4]),round.max=100);

                out.m[,5] <- round(as.numeric(out.m[,5]),2);

                write(paste(seedsSYM.v[m],
                            " (",nrow(seltopmodN.lm[[m]]),
                            " genes)",sep=""),
                      file=paste("topEpiModLists-",
                                 nameSTUDY,".txt",
                                 sep=""),
                      ncolumns=1,
                      append=TRUE);

                write.table(out.m,
                            file=paste("topEpiModLists-",
                                       nameSTUDY,
                                       ".txt",
                                       sep=""),
                            quote=FALSE,
                            sep="\t",
                            row.names=FALSE,
                            append=TRUE);
            }

        }

        return(list(size=sizeN.v,
                    mod=modN.v,
                    pv=modNpv.v,
                    selmod=selMod.idx,
                    fem=topmodN.m,
                    topmod=seltopmodN.lm,
                    sgc=sgcN.lo,
                    ew=tmpW.v,
                    adj=intEpi.o$adj));

    }

#' FemModShow
#'
#' https://rdrr.io/bioc/FEM/src/R/FemModShow.R
#'
#'
#' @param mod
#' @param name
#' @param fem.o
#' @param mode
#' @return A matrix of the infile
#' @export
FemModShow <-
    function(mod,
             name,
             fem.o,
             mode= "Integration"
             ){

        edgeweight=fem.o$ew
        adjacency=fem.o$adj

        ###############add shape to the vertex of igraph

        mycircle <- function(coords, v=NULL, params) {
            vertex.color <- params("vertex", "color")
            if (length(vertex.color) != 1 && !is.null(v)) {
                vertex.color <- vertex.color[v]
            }
            vertex.size  <- 1/200 * params("vertex", "size")
            if (length(vertex.size) != 1 && !is.null(v)) {
                vertex.size <- vertex.size[v]
            }
            vertex.frame.color <- params("vertex", "frame.color")
            if (length(vertex.frame.color) != 1 && !is.null(v)) {
                vertex.frame.color <- vertex.frame.color[v]
            }
            vertex.frame.width <- params("vertex", "frame.width")
            if (length(vertex.frame.width) != 1 && !is.null(v)) {
                vertex.frame.width <- vertex.frame.width[v]
            }

            mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
                   vertex.size, vertex.frame.width,
                   FUN=function(x, y, bg, fg, size, lwd) {
                       symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                               circles=size, add=TRUE, inches=FALSE)
                   })
        }

        #mod is from fem result such as HAND2,
        #hand2<-fembi.o$topmod$HAND2.
        #overall ideas: give the mtval, rtval,
        #vmcolor, vrcolor, label to vertex;
        #give weight, edgewidth to the edge of
        #realgraph. then subgraph
        #the HAND2 mod

        realgraph=graph.adjacency(adjacency,mode="undirected")

        #add the values
        E(realgraph)$weight= edgeweight #give the weight to the edges

        ############################################################
        # 1 edge width size

        edge.width.v=E(realgraph)$weight

        idxbw02=which(edge.width.v<=2 && edge.width.v>0)
        idxbw25=which(edge.width.v>2 && edge.width.v<5)
        idxbw510=which(edge.width.v>=5 && edge.width.v<10)
        idxgt10=which(edge.width.v>=10)

        edge.width.v[idxbw02]=1/4*edge.width.v[idxbw02]
        edge.width.v[idxbw25]=1/2*edge.width.v[idxbw25]
        edge.width.v[idxbw510]=3/4*edge.width.v[idxbw510]
        edge.width.v[idxgt10]=10#the max edgewidth is fixed as 10

        idxlt025=which(edge.width.v<0.25)
        edge.width.v[idxlt025]=0.25
        # if the edgewidth is less than 0.25,
        # it's too narrow to see. So fix the them with 0.25

        E(realgraph)$edgewidth=edge.width.v

        ##################################################################
        # 2 and 3  transform methylation value to node color, rna expression
        mod=as.data.frame(mod);
        if(mode=="Epi"){
            mod[,"stat(mRNA)"]=mod[,"stat(DNAm)"];
        }else if(mode=="Exp"){
            mod[,"stat(DNAm)"]=mod[,"stat(mRNA)"];
        }
        #if  there is mode is epi we add the stat(mRNA) colum and set them as 0;
        #sub graph mod and inhebited the edgewidth add the mval and rval
        mod.graph=igraph::induced.subgraph(realgraph,v=as.vector(mod[,1]))
        print(mod[,1])
        # the fembi.o$topmod$HAND2 is a dataframe and it has "Symbol",
        # "stat(DNAm)"
        # "stat(mRNA)"
        # which is useful
        mtval.v=vector()
        for(i in V(mod.graph)$name){
            mtval.v=c(mtval.v,(as.vector(mod[i,"stat(DNAm)"])))}
        # add the mtval to the mod graph one by one
        # according to the squence of mod graph name
        V(mod.graph)$mtval=mtval.v;

        rtval.v=vector()
        for(i in V(mod.graph)$name){
            rtval.v=c(rtval.v,(as.vector(mod[i,"stat(mRNA)"])))}
        # add the
        V(mod.graph)$rtval=rtval.v;
        print(rtval.v)

        # add the vm.color, vr.color
        vm.color=rep(0,length(V(mod.graph)));
        vr.color=rep(0,length(V(mod.graph)));
        # color scheme generate 100 colors
        tmcolor.scheme<-maPalette(low = "yellow",
                                  high ="blue",
                                  mid="grey",
                                  k =100);
        trcolor.scheme<-maPalette(low = "green",
                                  high ="red",
                                  mid="grey",
                                  k =100);

        tmcolor.scheme[13:88]="#BEBEBE";
        #( floor(-1.5/0.04)+51,floor(1.5/0.04)+51)
        #which is [13,88] is grey "#BEBEBE". 1.5 is
        #the thresh hold for the t values
        tmcolor.scheme[1:12]=maPalette(low = "yellow",
                                       high="lightyellow",
                                       k =12);
        tmcolor.scheme[89:100]=maPalette(low = "lightblue",
                                         high="blue",
                                         k =12);
        trcolor.scheme[13:88]="#BEBEBE";
        #[ floor(-1.5/0.04)+51,floor(1.5/0.04)+51]
        #which is [13,88] is grey "#BEBEBE". 1.5 is
        #the thresh hold for the t values
        trcolor.scheme[1:12]=maPalette(low = "green",
                                       high="lightgreen",
                                       k =12);
        trcolor.scheme[89:100]=maPalette(low = "#DC6868",
                                         high="red",
                                         k =12);
        # give the color according the mtval.
        # (-2,2),floor get the integer mod on 0.04 + 51
        # then get the according color.
        tmcolor.position=floor(as.numeric(V(mod.graph)$mtval)/0.04)+51;
        tmcolor.position[which(tmcolor.position<1)]<-1;
        tmcolor.position[which(tmcolor.position>100)]<-100;
        vm.color=tmcolor.scheme[tmcolor.position];
        V(mod.graph)$vmcolor<-vm.color ##add he vm.color to the vertex value
        print(vm.color)

        # add the frame color idea: get the tr color position then get
        # the color from trcolor.scheme
        trcolor.position=floor(as.numeric(V(mod.graph)$rtval)/0.04)+51;
        print(trcolor.position)
        trcolor.position[which(trcolor.position<1)]<-1;
        trcolor.position[which(trcolor.position>100)]<-100;
        vr.color=trcolor.scheme[trcolor.position];

        V(mod.graph)$vrcolor<-vr.color  ## the rna expression color

        if(mode=="Exp"){
            V(mod.graph)$color<-V(mod.graph)$vrcolor;
        }else{
            V(mod.graph)$color<-V(mod.graph)$vmcolor
            #use the vmcolor as the vertex color but if the mod is Exp
            #them the vertex color is vrcolor
        }
        print(vr.color)
        ####################################################################
        #add the mod label value

        label.v=vector()
        for(i in V(mod.graph)$name){
            label.v=c(label.v,(as.vector(mod[i,"Symbol"])))}
        #add the V(mod.graph)$name's
        #labels one by one from mod["$name","Symbol"]

        #####################################################################
        #create subgraph label.cex value
        V(mod.graph)$label.cex=rep(0.5,length(as.vector(V(mod.graph))));

        #all the cex first set as 0.7
        V(mod.graph)$label.cex[which(
            as.vector(V(mod.graph)$name)==as.vector(mod[1,1]))]=0.8
        #only the firt mod name was set as 1

        ####################################################################
        #generate the plot
        #when you want to plot the vertex shape,
        #and its frame width first you should load the api script api bellow
        add.vertex.shape("fcircle",
                         clip=igraph.shape.noclip,
                         plot=mycircle,
                         parameters=list(vertex.frame.color=1,
                                         vertex.frame.width=1))

        pdf(paste(name,".mod.pdf",sep=""))
        if(mode =="Integration"){

            plot(mod.graph,
                 layout=layout.fruchterman.reingold,
                 vertex.shape="fcircle",
                 vertex.frame.color=V(mod.graph)$vrcolor,
                 vertex.frame.width=4,
                 vertex.size=10,
                 vertex.label=label.v,
                 vertex.label.dist=0.6,
                 vertex.label.cex=V(mod.graph)$label.cex,
                 vertex.label.font=3,
                 edge.color="grey",
                 edge.width=E(mod.graph)$edgewidth)

            colorlegend(trcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(-0.5,0),
                        align="r",
                        cex=0.5)

            colorlegend(tmcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(0.5,1),
                        align="r",
                        cex=0.5)

            text(-1.50, 0.43, c("t(DNAm)\nCore"), cex=0.6)
            text(-1.50, -0.57,c("t(mRNA)\nBorder"), cex=0.6)

        }
        else if(mode =="Epi"){
            # if the mode is Epi the frame need not to show
            plot(mod.graph,
                 layout=layout.fruchterman.reingold,
                 vertex.frame.color=NA,
                 vertex.size=10,
                 vertex.label=label.v,
                 vertex.label.dist=0.6,
                 vertex.label.cex=V(mod.graph)$label.cex,
                 vertex.label.font=3,
                 edge.color="grey",
                 edge.width=E(mod.graph)$edgewidth)

            colorlegend(tmcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(0.5,1),
                        align="r",
                        cex=0.5)
            text(-1.50, 0.43, c("t(DNAm)"),cex=0.6)

        }
        else if(mode =="Exp"){

            plot(mod.graph,
                 layout=layout.fruchterman.reingold,
                 vertex.frame.color=NA,
                 vertex.size=10,
                 vertex.label=label.v,
                 vertex.label.dist=0.6,
                 vertex.label.cex=V(mod.graph)$label.cex,
                 vertex.label.font=3,
                 edge.color="grey",
                 edge.width=E(mod.graph)$edgewidth)

            colorlegend(trcolor.scheme,
                        seq(-2,2,0.5),
                        ratio.colbar=0.3,
                        xlim=c(-1.55,-1.4),
                        ylim=c(-0.5,0),
                        align="r",
                        cex=0.5)

            text(-1.50, -0.57,c("t(mRNA)"),cex=0.6)

        }
        dev.off()
        return(igraph.to.graphNEL(mod.graph));

}

#' champ.EpiMod
#'
#' BiocManager::install("FEM")
#' https://rdrr.io/github/gaberosser/ChAMP/src/R/champ.EpiMod.R
#'
#'
#' @param beta
#' @param pheno
#' @param nseeds
#' @param gamma
#' @param nMC
#' @param sizeR.v
#' @param minsizeOUT
#' @param resultsDir
#' @param PDFplot=TRUE
#' @param arraytype
#' @return A matrix of the infile
#' @export
champ.EpiMod <- function(beta=myNorm,
                         pheno=myLoad$pd$Sample_Group,
                         nseeds=100,
                         gamma=0.5,
                         nMC=1000,
                         sizeR.v=c(1,100),
                         minsizeOUT=10,
                         resultsDir="./CHAMP_EpiMod/",
                         PDFplot=TRUE,
                         arraytype="450K")
{

    if(getRversion() >= "3.1.0") utils::globalVariables(c("myNorm",
                                                        "myLoad",
                                                        "hprdAsigH.m")
    )
    message("[===========================]")
    message("[<<< ChAMP.EpiMod START >>>>]")
    message("-----------------------------")

    ### Prepare Checking ###
    if (!file.exists(resultsDir)) dir.create(resultsDir)
    message("champ.EpiMod Results will be saved in ",resultsDir)

    data(hprdAsigH)
    message("<< Load PPI network hprdAsigH >>")

    if(arraytype=="EPIC")
        statM.o <- GenStatM(beta,pheno,arraytype)
    else
        statM.o <- GenStatM(beta,pheno,"450k")
    message("<< Generate statM.o >>")

    intEpi.o=DoIntEpi450k(statM.o,hprdAsigH.m,c=1)

    message("<< Calculate EpiMod.o >>")
    EpiMod.o=DoEpiMod(intEpi.o,
                      nseeds=nseeds,
                      gamma=gamma,
                      nMC=nMC,
                      sizeR.v=sizeR.v,
                      minsizeOUT=minsizeOUT,
                      writeOUT=TRUE,
                      ew.v=NULL);

    if(PDFplot)
    {
        message("<< Draw All top significant module plot in PDF >>")
        tmpdir <- getwd()
        setwd(resultsDir)
        for(i in names(EpiMod.o$topmod)) FemModShow(EpiMod.o$topmod[[i]],
                                                    name=i,
                                                    EpiMod.o,
                                                    mode="Epi")
        setwd(tmpdir)
    }

    message("[<<<< ChAMP.EpiMod END >>>>>]")
    message("[===========================]")
    return(EpiMod.o)
}
