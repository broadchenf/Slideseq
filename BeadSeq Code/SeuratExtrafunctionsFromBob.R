

#A simple function to convert the weird 10x format into a sparse matrix for loading into Seurat.
Convert10xtoMatrix <-function(input){
  
  raw.data = exprs(input)
  genes = fData(input)
  idx.unique =which(table(genes[,2])==1)
  genes.use = names(table(genes[,2]))[idx.unique]
  idx = match(genes.use,genes[,2])
  raw.data = raw.data[idx,]
  rownames(raw.data) = genes.use
  return(raw.data)  
}

scalar1 <- function(x) {x / sqrt(sum(x^2))}


#I prefer how we select variable genes to Seurat's approach, which has some weird binning and normalization
mean.var.auto.seurat = function(object, alphathresh=0.99,varthresh=0.1,cex.use=0.3) {
  dge = object@raw.data
  xm = sweep(dge,2,colSums(dge),"/")
  trx_per_cell <- colSums(dge)
  
  gene_expr_mean <- rowMeans(xm)  # Each gene's mean expression level (across all cells)
  gene_expr_var  <- apply( xm, 1, var  )  # Each gene's expression variance (across all cells)
  nolan_constant <- mean((1 / trx_per_cell )) 
  alphathresh.corrected=alphathresh/dim(dge)[1]
  genemeanupper <- gene_expr_mean+qnorm(1-alphathresh.corrected/2)*sqrt(gene_expr_mean*nolan_constant/dim(dge)[2])
  genes.use=names(gene_expr_var)[which(gene_expr_var/nolan_constant> genemeanupper & log10(gene_expr_var) > log10(gene_expr_mean)+(log10(nolan_constant)+varthresh))]
  plot( log10(gene_expr_mean), log10(gene_expr_var), cex=cex.use
          )
  
  points(log10(gene_expr_mean[genes.use]),log10(gene_expr_var[genes.use]),cex=cex.use,col='green')
  abline(log10(nolan_constant),1,col='purple')
  
  legend("bottomright",paste0("Selected genes: ",length(genes.use)),pch=20,col='green')
  return(genes.use)
}


environment(mean.var.auto.seurat)<-asNamespace('Seurat')


#The main curation code for ICA.  Plots gene loadings and cell embeddings for each IC, and displays highest loading genes.  Helps determine whether an IC relates to a doublet class, an artifact, or is more likely to be a real biological signal.

#Intention is to run this twice for each ICA: run once with curation.data = NULL, which displays all ICs and draws a tSNE from all ICs; then, use the generated PDF to curate which ICs to retain, and feed in that curation data into the next function call.

CurateICA.seurat<-function(object,make.tsne = T, feature.plot=TRUE,kurtosis_threshold = 6,skewness_threshold = 1,curation.data=NULL,cluster.doublets = T,doubletalpha=0.01) {
  
  x = list(object@dr$ica@gene.loadings,object@dr$ica@cell.embeddings)
  
  n_ICs = dim(x[[1]])[2]
  # Calculate skews of cell scores and gene loadings
  gene_skews_initial = apply( x[[1]], 2, skewness )
  cell_skews_initial = apply( x[[2]], 2, skewness )
  
  # Flip certain ICs so that all ICs have positive skew for cell scores
  wh_ICs_flip = which( cell_skews_initial < 0 )
  x[[1]][,wh_ICs_flip] = -1 * x[[1]][,wh_ICs_flip]
  x[[2]][,wh_ICs_flip] = -1 * x[[2]][,wh_ICs_flip]
  
  # Reorder the ICs by their skew (from best to worst)
  IC_reorder = order( abs(cell_skews_initial), decreasing=T )
  x[[1]] = x[[1]][,IC_reorder]
  x[[2]] = x[[2]][,IC_reorder]
  
  # Calculate cell-kurtosis(k), gene-skew(gs), and cell-skew(cs) for the reordered ICs
  gk = apply( x[[1]], 2, kurtosis )
  ck = apply( x[[2]], 2, kurtosis )
  gs = apply( x[[1]], 2, skewness )
  cs = apply( x[[2]], 2, skewness )
  if (!is.null(curation.data)) {
  curation.names = curation.data$status
  
  keep_IC=(curation.names == "Real")
  
  idx.keep=which(keep_IC)
  
  } else {
    keep_IC = rep(TRUE,times = ncol(x[[2]]))
    idx.keep = which(keep_IC)
    cluster.doublets = F
    curation.names = rep("Uncurated",times = ncol(x[[2]]))
  }
  
  
  
  # Routine for 1-D clustering to remove Doublet and Outlier cells.
  if(cluster.doublets) {
    
    
    idx.doublet=which(grepl(curation.names,pattern="Doublet",ignore.case=T) | grepl(curation.names,pattern="Outlier",ignore.case=T))
    
    #Obtain the mode for centering the 1-D gaussian using the density() function, and the turnpoints function from pastecs
    
    doublets=lapply(idx.doublet,function(i){
      q=x[[2]][,i]
      d = density(q)
      ts_y = ts(d$y)
      tp = turnpoints(ts_y)
      idx = which.max(d$y[tp$tppos])
      center=d$x[tp$tppos[idx]]
      pvalsUpper=pnorm(q,mean=center,sd= sd(q),lower.tail=F)
      pvalsUpper=p.adjust(pvalsUpper,"fdr")
      idx=which(pvalsUpper <= doubletalpha)
      return(idx)
    })
    
    
  } else {
    idx.doublet=integer(0)
  }

  #Store reoriented ICs
  colnames(x[[2]])<-paste0("IC",1:ncol(x[[2]]))
  colnames(x[[1]])<-paste0("IC",1:ncol(x[[2]]))
  
  
  object@dr$ica@cell.embeddings = x[[2]]
  object@dr$ica@gene.loadings = x[[1]]
  object@dr$ica@misc = curation.names
  
  #If plotting on new tSNE, make new embedding with only retained ICs.
  if (make.tsne)
  object = RunTSNE(object = object,reduction.use = 'ica',dims.use = which(keep_IC),do.fast=T,pca=F,check_duplicates=F)

  
  # Skew-kurtosis plot of the ICs
  par(mfrow=c(1,1))
  plot(  log(ck,2), cs , col="white",
         xlab = "log2( Cell-score kurtosis)", ylab="Cell-score skew",
         main = "Cell-score skew and kurtosis (after component reordering)") 
  rect( 0, 0, log2(kurtosis_threshold), skewness_threshold, col="lightgray", border=NA )


  text( log(ck, 2), cs, 1:length(cs)  , col="black" )
  
  abline(v=0); abline(h=0)
  
  #  Review the cell-score and gene-loading distributions for each component

  for(i in 1:n_ICs) {
    par(mfrow=c(2,1))
    print(i)
    top_genes = row.names( x[[1]] )[ order(x[[1]][,i], decreasing=T )[1:7] ]
    top_genes_textstring = paste(top_genes, collapse=", ")
    stats_textstring1 = paste( "k=", round(ck[i],1), ", s=", round(cs[i],1))
    stats_textstring2 = paste( "k=", round(gk[i],1), ", s=", round(gs[i],1))
    component_textstring = paste( "Component", i, "  (", curation.names[i],")" )
    if( !keep_IC[i] ) component_textstring=paste(component_textstring, "(REMOVE)" )
    plot_title1 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring1 )
    plot_title2 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring2 )
    
    # Plot cell scores
    plot( x[[2]][,i], cex=0.2, main=plot_title1, xlab="Cell", ylab="Score"  )
    
    # Highlight the doublet-flagged cells.
    if (i %in% idx.doublet) {
      idx.use = which(idx.doublet==i)

      points(doublets[[idx.use]],x[[2]][doublets[[idx.use]],i],cex=0.5,col='red')
    }
    
    # Plot gene loadings
    plot( x[[1]][,i], cex=0.2, main=plot_title2, xlab="Gene", ylab="Loading"  )
    if (feature.plot) { 
      tsrots = object@dr$tsne@cell.embeddings
      par(mfrow=c(1,1))
      fplot(tsrots[rownames(x[[2]]),],x[[2]][,i],title=component_textstring)
      }
  }
  
  # Skew-skew plot of the ICs
  par(mfrow=c(1,1))
  plot(  gs, cs , col="white",
         xlab = "Gene-loading skew", ylab="Cell-score skew",
         main = "Cell-score and gene-loading skew (after component reordering)") 
  text( gs, cs,1:length(cs)   , col="black" )
  abline(v=0); abline(h=0)
  

  if(cluster.doublets) {
    object = SubsetData(object = object,cells.use = setdiff(rownames(x[[2]]),rownames(x[[2]])[unlist(doublets)]),do.scale = F)
    object@raw.data=object@raw.data[,colnames(object@data)]
    object@raw.data=object@raw.data[which(Matrix::rowSums(object@raw.data)>0),]
  }
    return(object)
}

fplot<-function(tsne,factor,title,cols.use=heat.colors(10),pt.size=0.7,pch.use=20) {
  c=intersect(rownames(tsne),names(factor))
  data.use=factor[c]
  data.cut=as.numeric(as.factor(cut(as.numeric(data.use),breaks=length(cols.use))))
  data.col=rev(cols.use)[data.cut]
  plot(tsne[c,1],tsne[c,2],col=data.col,cex=pt.size,pch=pch.use,main=title)
  
}

environment(CurateICA.seurat)<-asNamespace('Seurat')

StoreSpatialData<-function(object,spatial_coord) {
  cells.use = intersect(names(object@ident),rownames(spatial_coord))
  spatial.obj = new(Class="dim.reduction",gene.loadings=matrix(0),cell.embeddings=spatial_coord[cells.use,],key="SPATIAL_")
  colnames(spatial.obj@cell.embeddings) = paste0("SPATIAL",1:ncol(spatial.obj@cell.embeddings)) 
  object = SubsetData(object,cells.use = cells.use)
  object@raw.data = object@raw.data[,colnames(object@data)]
  object@dr$spatial=spatial.obj
  return(object)
}
environment(StoreSpatialData)<-asNamespace('Seurat')



AtlasToSeurat<-function(data_folder_path,subcluster=NULL) {
    data_folder_path=sub(x=data_folder_path,pattern = "/$",replacement = "")
    
    split=strsplit(data_folder_path,split="/")[[1]]
    samplename=split[length(split)]
    idents=readRDS(paste0(data_folder_path,"/assign/",samplename,".subcluster.assign.RDS"))
    if (is.null(subcluster)) {
        rawdata=readRDS(paste0(data_folder_path,"/dge/",samplename,".filtered.raw.dge.RDS"))
        normalized.dge=readRDS(paste0(data_folder_path,"/dge/",samplename,".filtered.scaled.dge.RDS"))
        orig.dge=rawdata[,names(idents)]
        tsne=readRDS(paste0(data_folder_path,"/tSNE/",samplename,"_tSNExy.RDS"))
        ics=readRDS(paste0(data_folder_path,"/components/",samplename,".ica.RDS"))
        my.seurat=CreateSeuratObject(raw.data=orig.dge)
        my.seurat@data=normalized.dge
     
        my.seurat@var.genes = rownames(ics$gene_loadings)
        my.seurat@scale.data=t(scale(t(normalized.dge[my.seurat@var.genes,]),center=T,scale=T))
        ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations,key="IC")
        colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))    
        colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
        tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
        colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
        cells.int=intersect(rownames(tsne),names(idents))
        cells.removed=setdiff(rownames(tsne),cells.int)
        idents.new=factor(x = c(as.character(idents[cells.int]),rep("REMOVED",times=length(cells.removed))))
        names(idents.new) = c(cells.int,cells.removed)
        idents.new = idents.new[rownames(tsne)]
        my.seurat@dr$ica=ica.obj
        my.seurat@dr$tsne=tsne.obj
        my.seurat@ident=idents.new
        return(my.seurat)
        
        
        
    } else {
        rawdata=readRDS(paste0(data_folder_path,"/cluster",subcluster,"/",samplename,".subcluster_inputs.RDS"))
        orig.dge=rawdata$raw_dge
        scaled.dge=rawdata$scaled_gene_selected_dge
        
        tsne=readRDS(paste0(data_folder_path,"/tSNE/",samplename,".cluster",subcluster,".CURATEDtSNE.RDS"))
        if(is.null(tsne)) {
            tsne=readRDS(paste0(data_folder_path,"/cluster",subcluster,"/",samplename,".cluster",subcluster,".auto.tSNExy.RDS"))
            
        }
        icfiles=list.files(paste0(data_folder_path,"/components/"))
        icfile.use=icfiles[grep(icfiles,pattern=paste0("cluster",subcluster,"\\."))]
        ics=readRDS(paste0(data_folder_path,"/components/",icfile.use))
        curationsheets=list.files(paste0(data_folder_path,"/curation_sheets/"))
        curationfile=curationsheets[grep(curationsheets,pattern=paste0("Subcluster_",subcluster,"[_\\.]"))]
        ic.annotations = read.table(paste0(data_folder_path,"/curation_sheets/",curationfile), stringsAsFactors=F, sep=",", header=T, quote="\"")
        my.seurat=CreateSeuratObject(raw.data=orig.dge)
     
        my.seurat@data = sweep(my.seurat@raw.data,2,Matrix::colSums(my.seurat@raw.data),"/")
        my.seurat@scale.data=scaled.dge
        my.seurat@var.genes = rownames(ics$gene_loadings)
        ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations,key="IC",misc=ic.annotations)
        colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))
        colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
        tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
        colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
        cells.int=intersect(rownames(tsne),names(idents))
        cells.removed=setdiff(rownames(tsne),cells.int)
        idents.new=factor(x = c(as.character(idents[cells.int]),rep("REMOVED",times=length(cells.removed))))
        names(idents.new) = c(cells.int,cells.removed)
        idents.new = idents.new[rownames(tsne)]
        my.seurat@dr$ica=ica.obj
        my.seurat@dr$tsne=tsne.obj
        my.seurat@ident=idents.new
        return(my.seurat)
        
    }
    
}

AtlasToSeuratWindows<-function(data_folder_path,subcluster=NULL) {
  data_folder_path=sub(x=data_folder_path,pattern = "\\$",replacement = "")
  
  split=strsplit(data_folder_path,split="\\\\")[[1]]
  samplename=split[length(split)]
  idents=readRDS(paste0(data_folder_path,"\\assign\\",samplename,".subcluster.assign.RDS"))
  if (is.null(subcluster)) {
    rawdata=readRDS(paste0(data_folder_path,"\\dge\\",samplename,".filtered.raw.dge.RDS"))
    normalized.dge=readRDS(paste0(data_folder_path,"\\dge\\",samplename,".filtered.scaled.dge.RDS"))
    orig.dge=rawdata[,names(idents)]
    tsne=readRDS(paste0(data_folder_path,"\\tSNE\\",samplename,"_tSNExy.RDS"))
    ics=readRDS(paste0(data_folder_path,"\\components\\",samplename,".ica.RDS"))
    my.seurat=CreateSeuratObject(raw.data=orig.dge)
    my.seurat@data=normalized.dge
    
    my.seurat@var.genes = rownames(ics$gene_loadings)
    my.seurat@scale.data=t(scale(t(normalized.dge[my.seurat@var.genes,]),center=T,scale=T))
    ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations,key="IC")
    colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))    
    colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
    tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
    cells.int=intersect(rownames(tsne),names(idents))
    cells.removed=setdiff(rownames(tsne),cells.int)
    idents.new=factor(x = c(as.character(idents[cells.int]),rep("REMOVED",times=length(cells.removed))))
    names(idents.new) = c(cells.int,cells.removed)
    idents.new = idents.new[rownames(tsne)]
    my.seurat@dr$ica=ica.obj
    my.seurat@dr$tsne=tsne.obj
    my.seurat@ident=idents.new
    return(my.seurat)
    
    
    
  } else {
    rawdata=readRDS(paste0(data_folder_path,"\\cluster",subcluster,"\\",samplename,".subcluster_inputs.RDS"))
    orig.dge=rawdata$raw_dge
    scaled.dge=rawdata$scaled_gene_selected_dge
    
    tsne=readRDS(paste0(data_folder_path,"\\tSNE\\",samplename,".cluster",subcluster,".CURATEDtSNE.RDS"))
    if(is.null(tsne)) {
      tsne=readRDS(paste0(data_folder_path,"\\cluster",subcluster,"\\",samplename,".cluster",subcluster,".auto.tSNExy.RDS"))
      
    }
    icfiles=list.files(paste0(data_folder_path,"\\components\\"))
    icfile.use=icfiles[grep(icfiles,pattern=paste0("cluster",subcluster,"\\."))]
    ics=readRDS(paste0(data_folder_path,"\\components\\",icfile.use))
    curationsheets=list.files(paste0(data_folder_path,"\\curation_sheets\\"))
    curationfile=curationsheets[grep(curationsheets,pattern=paste0("Subcluster_",subcluster,"[_\\.]"))]
    ic.annotations = read.table(paste0(data_folder_path,"\\curation_sheets\\",curationfile), stringsAsFactors=F, sep=",", header=T, quote="\"")
    my.seurat=CreateSeuratObject(raw.data=orig.dge)
    
    my.seurat@data = sweep(my.seurat@raw.data,2,Matrix::colSums(my.seurat@raw.data),"\\")
    my.seurat@scale.data=scaled.dge
    my.seurat@var.genes = rownames(ics$gene_loadings)
    ica.obj = new(Class="dim.reduction",gene.loadings=ics$gene_loadings,cell.embeddings=ics$cell_rotations,key="IC",misc=ic.annotations)
    colnames(ica.obj@cell.embeddings) = paste0("IC",1:ncol(ica.obj@cell.embeddings))
    colnames(ica.obj@gene.loadings) = paste0("IC",1:ncol(ica.obj@gene.loadings))
    tsne.obj=new(Class="dim.reduction",cell.embeddings=tsne,key="tSNE_")
    colnames(tsne.obj@cell.embeddings) = paste0("tSNE_",1:ncol(tsne.obj@cell.embeddings))
    cells.int=intersect(rownames(tsne),names(idents))
    cells.removed=setdiff(rownames(tsne),cells.int)
    idents.new=factor(x = c(as.character(idents[cells.int]),rep("REMOVED",times=length(cells.removed))))
    names(idents.new) = c(cells.int,cells.removed)
    idents.new = idents.new[rownames(tsne)]
    my.seurat@dr$ica=ica.obj
    my.seurat@dr$tsne=tsne.obj
    my.seurat@ident=idents.new
    return(my.seurat)
    
  }
  
}


RunNMF<-function(object, factors.compute = 50, print.results = TRUE, factors.print = 1:factors.compute, genes.print = 50, seed.use = 1, ...) {
     norm=log(sweep(object@raw.data,2,colSums(object@raw.data),"/")*10000+1)
     #norm=sweep(object@raw.data,2,colSums(object@raw.data),"/")
     
    uncentered.scaled = scale(t(norm[object@var.genes,]),center=F,scale=T)
     # uncentered.scaled = apply(norm[object@var.genes,],1,scalar1)
    factors.compute=min(c(factors.compute,dim(uncentered.scaled)))
    set.seed(seed=seed.use)
    nmf.results = nnmf(uncentered.scaled,k=factors.compute)
    nmf.obj<-new(Class="dim.reduction",gene.loadings=t(nmf.results$H),cell.embeddings=nmf.results$W,key="NMF")
    colnames(nmf.obj@cell.embeddings) = paste0("Factor",1:ncol(nmf.obj@cell.embeddings))
    colnames(nmf.obj@gene.loadings) = paste0("Factor",1:ncol(nmf.obj@gene.loadings))
    object@dr$nmf <- nmf.obj
    if(print.results) {
        for (i in factors.print) {
            code <- paste0(GetDimReduction(object = object, reduction.type = "nmf",
            slot = "key"), i)
            sx <- DimTopGenes(object = object, dim.use = i, reduction.type = 'nmf',
            num.genes = genes.print*2, use.full = F,do.balanced = FALSE)
            print(code)
            print((sx[1:genes.print]))
            print("")
            print(rev(x = (sx[(length(x = sx) - genes.print + 1):length(x = sx)])))
            print("")
            print("")
        }
        
        
    }
    return(object)
}

environment(RunNMF)<-asNamespace('Seurat')



CurateNMF.seurat<-function(object,make.tsne = T, feature.plot=TRUE,kurtosis_threshold = 6,skewness_threshold = 1,curation.data=NULL,cluster.doublets = T,doubletalpha=0.01) {
  
  x = list(object@dr$nmf@gene.loadings,object@dr$nmf@cell.embeddings)
  
  n_ICs = dim(x[[1]])[2]
  # Calculate skews of cell scores and gene loadings
  gene_skews_initial = apply( x[[1]], 2, skewness )
  cell_skews_initial = apply( x[[2]], 2, skewness )
  
  # Flip certain ICs so that all ICs have positive skew for cell scores
  wh_ICs_flip = which( cell_skews_initial < 0 )
  x[[1]][,wh_ICs_flip] = -1 * x[[1]][,wh_ICs_flip]
  x[[2]][,wh_ICs_flip] = -1 * x[[2]][,wh_ICs_flip]
  
  # Reorder the ICs by their skew (from best to worst)
  IC_reorder = order( abs(cell_skews_initial), decreasing=T )
  x[[1]] = x[[1]][,IC_reorder]
  x[[2]] = x[[2]][,IC_reorder]
  
  # Calculate cell-kurtosis(k), gene-skew(gs), and cell-skew(cs) for the reordered ICs
  gk = apply( x[[1]], 2, kurtosis )
  ck = apply( x[[2]], 2, kurtosis )
  gs = apply( x[[1]], 2, skewness )
  cs = apply( x[[2]], 2, skewness )
  if (!is.null(curation.data)) {
    curation.names = curation.data$status
    
    keep_IC=(curation.names == "Real")
    
    idx.keep=which(keep_IC)
    
  } else {
    keep_IC = rep(TRUE,times = ncol(x[[2]]))
    idx.keep = which(keep_IC)
    cluster.doublets = F
    curation.names = rep("Uncurated",times = ncol(x[[2]]))
  }
  
  
  
  # Routine for 1-D clustering to remove Doublet and Outlier cells.
  if(cluster.doublets) {
    
    
    idx.doublet=which(grepl(curation.names,pattern="Doublet",ignore.case=T) | grepl(curation.names,pattern="Outlier",ignore.case=T))
    
    #Obtain the mode for centering the 1-D gaussian using the density() function, and the turnpoints function from pastecs
    
    doublets=lapply(idx.doublet,function(i){
      q=x[[2]][,i]
      d = density(q)
      ts_y = ts(d$y)
      tp = turnpoints(ts_y)
      idx = which.max(d$y[tp$tppos])
      center=d$x[tp$tppos[idx]]
      pvalsUpper=pnorm(q,mean=center,sd= sd(q),lower.tail=F)
      pvalsUpper=p.adjust(pvalsUpper,"fdr")
      idx=which(pvalsUpper <= doubletalpha)
      return(idx)
    })
    
    
  } else {
    idx.doublet=integer(0)
  }
  
  #Store reoriented ICs
  colnames(x[[2]])<-paste0("Factor",1:ncol(x[[2]]))
  colnames(x[[1]])<-paste0("Factor",1:ncol(x[[2]]))
  object@dr$nmf@cell.embeddings = x[[2]]
  object@dr$nmf@gene.loadings = x[[1]]
  object@dr$nmf@misc = curation.names
  
  #If plotting on new tSNE, make new embedding with only retained ICs.
  if (make.tsne)
    object = RunTSNE(object = object,reduction.use = 'ica',dims.use = which(keep_IC),do.fast=T,pca=F)
  
  
  # Skew-kurtosis plot of the ICs
  par(mfrow=c(1,1))
  plot(  log(ck,2), cs , col="white",
         xlab = "log2( Cell-score kurtosis)", ylab="Cell-score skew",
         main = "Cell-score skew and kurtosis (after component reordering)") 
  rect( 0, 0, log2(kurtosis_threshold), skewness_threshold, col="lightgray", border=NA )
  
  
  text( log(ck, 2), cs, 1:length(cs)  , col="black" )
  
  abline(v=0); abline(h=0)
  
  #  Review the cell-score and gene-loading distributions for each component
  
  for(i in 1:n_ICs) {
    par(mfrow=c(2,1))
    print(i)
    top_genes = row.names( x[[1]] )[ order(x[[1]][,i], decreasing=T )[1:7] ]
    top_genes_textstring = paste(top_genes, collapse=", ")
    stats_textstring1 = paste( "k=", round(ck[i],1), ", s=", round(cs[i],1))
    stats_textstring2 = paste( "k=", round(gk[i],1), ", s=", round(gs[i],1))
    component_textstring = paste( "Component", i, "  (", curation.names[i],")" )
    if( !keep_IC[i] ) component_textstring=paste(component_textstring, "(REMOVE)" )
    plot_title1 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring1 )
    plot_title2 = paste( component_textstring, "\n", top_genes_textstring, "\n", stats_textstring2 )
    
    # Plot cell scores
    plot( x[[2]][,i], cex=0.2, main=plot_title1, xlab="Cell", ylab="Score"  )
    
    # Highlight the doublet-flagged cells.
    if (i %in% idx.doublet) {
      idx.use = which(idx.doublet==i)
      
      points(doublets[[idx.use]],x[[2]][doublets[[idx.use]],i],cex=0.5,col='red')
    }
    
    # Plot gene loadings
    plot( x[[1]][,i], cex=0.2, main=plot_title2, xlab="Gene", ylab="Loading"  )
    if (feature.plot) { 
      tsrots = object@dr$tsne@cell.embeddings
      par(mfrow=c(1,1))
      fplot(tsrots[rownames(x[[2]]),],x[[2]][,i],title=component_textstring)
    }
  }
  
  # Skew-skew plot of the ICs
  par(mfrow=c(1,1))
  plot(  gs, cs , col="white",
         xlab = "Gene-loading skew", ylab="Cell-score skew",
         main = "Cell-score and gene-loading skew (after component reordering)") 
  text( gs, cs,1:length(cs)   , col="black" )
  abline(v=0); abline(h=0)
  
  
  if(cluster.doublets) {
    object = SubsetData(object = object,cells.use = setdiff(rownames(x[[2]]),unlist(doublets)),do.scale = F)
  }
  return(object)
}
environment(CurateNMF.seurat)<-asNamespace('Seurat')


fplot<-function(tsne,factor,title,cols.use=heat.colors(10),pt.size=0.7,pch.use=20) {
  c=intersect(rownames(tsne),names(factor))
  data.use=factor[c]
  data.cut=as.numeric(as.factor(cut(as.numeric(data.use),breaks=length(cols.use))))
  data.col=rev(cols.use)[data.cut]
  plot(tsne[c,1],tsne[c,2],col=data.col,cex=pt.size,pch=pch.use,main=title)
  
}

correlate.factors<-function(object1,object2) {
  H1=object1@dr$nmf@gene.loadings
  H1=H1[,which(object1@dr$nmf@misc=="Real")]
  H2=object2@dr$nmf@gene.loadings
  H2=H2[,which(object2@dr$nmf@misc=="Real")]
  shared.genes = intersect(rownames(H1),rownames(H2))
  cors=cor(H1[shared.genes,],H2[shared.genes,])
  if(ncol(H1)==ncol(H2))
  cors=matrix.sort(cors)
  return(cors)
  
}

correlate.ICs<-function(object1,object2) {
  H1=object1@dr$ica@gene.loadings
  H2=object2@dr$ica@gene.loadings
  shared.genes = intersect(rownames(H1),rownames(H2))
  cors=cor(H1[shared.genes,],H2[shared.genes,])
  if(ncol(H1)==ncol(H2))
    cors=matrix.sort(cors)
  return(cors)
  
}



matrix.sort <- function(matrix, require_square=TRUE) {
  
  if (require_square && nrow(matrix) != ncol(matrix)) stop("Not diagonal")
  if(is.null(rownames(matrix))) rownames(matrix) <- 1:nrow(matrix)
  
  row.max <- apply(matrix,1,which.max)
  if(all(table(row.max) != 1)) stop("Ties cannot be resolved")
  
  matrix[names(sort(row.max)),]
}


BuildRFClassifierSCALED<- function (object, training.genes = NULL, training.classes = NULL, 
          verbose = TRUE, ...) 
{
  training.classes <- as.vector(x = training.classes)
  training.genes <- SetIfNull(x = training.genes, default = rownames(x = object@data))
  training.data <- as.data.frame(x = as.matrix(x = scale(t(x = object@data[training.genes, 
                                                                     ])),center=T,scale=T))
  training.data$class <- factor(x = training.classes)
  if (verbose) {
    print("Training Classifier ...")
  }
  classifier <- ranger(data = training.data, dependent.variable.name = "class", 
                       classification = TRUE, write.forest = TRUE, ...)
  return(classifier)
}
environment(BuildRFClassifierSCALED)<-asNamespace('Seurat')


downsample.matrix <- function(mat, samplerate=0.8,seed=1) {
    set.seed(seed)
    new <- matrix(0, nrow(mat), ncol(mat))
    colnames(new) <- colnames(mat)
    rownames(new) <- rownames(mat)
    for (i in 1:nrow(mat)) {
        for (j in 1:ncol(mat)) {
            new[i,j] <- sum(runif(mat[i,j], 0, 1) < samplerate)
        }
    }
    return(new)
}


fastRead<-function(inFile) {
  if (length(grep(".gz", inFile))==1) {
    t=tempfile()
    cmd=paste("gunzip -c", inFile, ">", t)  
    system(cmd,intern=F)	
    inFileFinal=t
  } else {
    inFileFinal=inFile
  }
  
  
  a=fread(inFileFinal, data.table=F,skip = 3,header = T)
  rownames(a)<-a$GENE
  a=a[,-1]
  return (a)
}

