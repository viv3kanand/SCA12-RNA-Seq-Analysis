# Read salon quant files using tximport for DGE and DTE
tximport = function(dir, coldata, txdf, DTE=TRUE){
  require(tximport)
  require(SummarizedExperiment)
  
  files <- file.path(dir, list.files(dir), "quant.sf")
  names(files) <- list.files(dir)
  
  if(all(file.exists(files))){ # check if all files exist
    print("All files exist in the salmon folder")
  }
  
  if(DTE){
    txi <- tximport(files, type = "salmon", tx2gene = txdf, txOut = TRUE)
  }
  else {
    txi <- tximport(files, type = "salmon", tx2gene = txdf, countsFromAbundance = "lengthScaledTPM")
    # alternative usage -> tximport::summarizeToGene(txi, txdf, countsFromAbundance = "lengthScaledTPM")
  }
  
  coldata <- merge(x = data.frame(names = names(files),
                                  files = files), 
                   y = coldata,
                   by.x = "names", by.y = "Sample", 
                   all=TRUE, sort = FALSE)
  
  tximportInfo <- list(version=packageVersion("tximport"),
                      type=type,
                      importTime=Sys.time())
  metadata <- list(tximportInfo=tximportInfo)
  metadata$countsFromAbundance <- txi$countsFromAbundance
  
  se <- makeUnrangedSE(txi, coldata, metadata, txdf)
  
  return(se)
}
  
  
#### from tximeta wrapper for adding Inferential replicates in se because of dplyr filter error ########
makeUnrangedSE <- function(txi, coldata, metadata, txdf) {
  assays <- txi[c("counts","abundance","length")]
  # if there are inferential replicates
  if ("infReps" %in% names(txi)) {
    infReps <- rearrangeInfReps(txi$infReps)
    assays <- c(assays, infReps)
  }
  
  SummarizedExperiment(assays=assays,
                       colData=coldata,
                       metadata=metadata)
}


rearrangeInfReps <- function(infReps) {
  nreps <- ncol(infReps[[1]])
  stopifnot(all(sapply(infReps, ncol) == nreps))
  getCols <- function(j,l) do.call(cbind, lapply(seq_along(l), function(k)  l[[k]][,j]))
  infReps <- lapply(seq_len(nreps), getCols, infReps)
  names(infReps) <- paste0("infRep",seq_len(nreps))
  infReps
}
######################################################################################################

# Replicating DESeq2::plotPCA top variable genes for df and vst
call_PCA <- function(x, ntop=500){
  require(SummarizedExperiment)
  y <- if(is.data.frame(x)) as.matrix(x) else assay(x)
  rv <- rowVars(y)
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  
  pca <- prcomp(t(y[select,]), center = TRUE, scale. = FALSE)
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  return(list("pca" = pca, "percentVar" = percentVar))
}























