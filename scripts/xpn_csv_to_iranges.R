#!/usr/bin/Rscript

qw <- function(...) {
  as.character(sys.call()[-1])
}

options(stringsAsFactors = FALSE);

# load datafile from pp_expression_data
# In this case, we have the F test results, and the
# individual t-test contrasts
# 1: BMP v CTX
# 2: FBS v CTX
# 3: BMP v FBS

f.test <- read.csv('/mnt/data/limma_results.csv')
t.test.1 <- read.csv('/mnt/data/limma_results_bmp_v_ctx.csv')
t.test.2 <- read.csv('/mnt/data/limma_results_fbs_v_ctx.csv')
t.test.3 <- read.csv('/mnt/data/limma_results_bmp_v_fbs.csv')

f.test <- f.test[,-1]
t.test.1 <- t.test.1[,-1]
t.test.2 <- t.test.2[,-1]
t.test.3 <- t.test.3[,-1]
colnames(f.test)[1] <- colnames(t.test.1)[1] <- colnames(t.test.2)[1] <- colnames(t.test.3)[1] <- "IlluminaID"

remoat.annot.file <- "/mnt/work/lib/Annotation_Illumina_Mouse-WG-V2_mm9_V1.0.0_Aug09.txt"




# Array is Illumina Mouse Ref8 v2
# "This BeadChip targets approximately 25,600 well-annotated RefSeq transcripts"
#dim(limma)
#[1] 25697     7


# The IDs we have are bead IDs not ILMN ids. Just seems to be what GIS give us back.
# We can annotate these using
library(illuminaMousev2BeadID.db)
# but this library doesn't have the genome position of the probes. 
# The BioC annotation comes from the ReMOAT group at Cambridge Uni, so
# we can just directly use their data.
# ReMOAT paper: http://nar.oxfordjournals.org/cgi/content/full/gkp942

remoat.annot <- read.csv(remoat.annot.file, header=T, sep="\t") 

# Throw away anything that isn't perfect or good (the only classes of match likely
# to give a reliable signal according to the ReMOAT paper.
classes <- unique(remoat.annot[,"Quality_score"])
perfect <- grep("Perfect", classes)
good <- grep("Good", classes)
classes <- classes[union(perfect, good)]
keep <- which(remoat.annot[,"Quality_score"] %in% classes)
remoat.annot <- remoat.annot[keep,]


# join illumina data to annotation. This will just ditch any limma data from unmapped probes
# (limma = 15060, limma.annot=13899)
f.test.annot <- merge(f.test, remoat.annot, by.x="IlluminaID", by.y="Array_Address_Id_0")
t.test.1.annot <- merge(t.test.1, remoat.annot, by.x="IlluminaID", by.y="Array_Address_Id_0")
t.test.2.annot <- merge(t.test.2, remoat.annot, by.x="IlluminaID", by.y="Array_Address_Id_0")
t.test.3.annot <- merge(t.test.3, remoat.annot, by.x="IlluminaID", by.y="Array_Address_Id_0")

# We probably don't need *all* the annotation, although we do need to have the genome and
# transcriptome position of the probe.
# Multiple definitions of transcriptome though, so we need results for all of them.
# see http://www.compbio.group.cam.ac.uk/Resources/Annotation/Description_1.0.0.xls
# for colname definitions


f.nms <- qw(IlluminaID, bmpvctx, fbsvctx, bmpvfbs, AveExpr, F,P.Value, adj.P.Val,
          Lumi_id,Probe_sequence,Probe_type, Quality_score, Genomic_location,
          RefSeq_transcripts, Proportion_RefSeq_transcripts,
          UCSC_transcripts,Proportion_UCSC_transcripts,
          GenBank_transcripts, Proportion_GenBank_transcripts,
          Exons,
          Ensembl_transcripts, Proportion_Ensembl_transcripts,
          Lumi_transcriptomic_annotation, Lumi_transcriptomic_match )
f.test.annot <- f.test.annot[,f.nms]

t.nms <- qw(IlluminaID, logFC, AveExpr, t, P.Value, adj.P.Val, B,Lumi_id,Probe_sequence,Probe_type, Quality_score, Genomic_location,
          RefSeq_transcripts, Proportion_RefSeq_transcripts,
          UCSC_transcripts,Proportion_UCSC_transcripts,
          GenBank_transcripts, Proportion_GenBank_transcripts,
          Exons,
          Ensembl_transcripts, Proportion_Ensembl_transcripts,
          Lumi_transcriptomic_annotation, Lumi_transcriptomic_match )

t.test.1.annot <- t.test.1.annot[,t.nms]
t.test.2.annot <- t.test.2.annot[,t.nms]
t.test.3.annot <- t.test.3.annot[,t.nms]

#and while we're at it, we might at well get the Bioconductor annotation
ids <- as.character(f.test.annot[,"IlluminaID"])
entrez <-  unlist(mget(ids,illuminaMousev2BeadIDENTREZID, ifnotfound=NA))
ensembl <- unlist(mget(ids,illuminaMousev2BeadIDENSEMBL, ifnotfound=NA ))
genename <- unlist(mget(ids,illuminaMousev2BeadIDGENENAME, ifnotfound=NA))
genesymbol <- unlist(mget(ids,illuminaMousev2BeadIDSYMBOL, ifnotfound=NA))

f.test.annot <- data.frame(f.test.annot, EntrezGene=entrez[ids], EnsemblGene=ensembl[ids], GeneName=genename[ids], GeneSymbol=genesymbol[ids])

ids <- as.character(t.test.1.annot[,"IlluminaID"])
t.test.1.annot <- data.frame(t.test.1.annot, EntrezGene=entrez[ids], EnsemblGene=ensembl[ids], GeneName=genename[ids], GeneSymbol=genesymbol[ids])

ids <- as.character(t.test.2.annot[,"IlluminaID"])
t.test.2.annot <- data.frame(t.test.2.annot, EntrezGene=entrez[ids], EnsemblGene=ensembl[ids], GeneName=genename[ids], GeneSymbol=genesymbol[ids])

ids <- as.character(t.test.3.annot[,"IlluminaID"])
t.test.3.annot <- data.frame(t.test.3.annot, EntrezGene=entrez[ids], EnsemblGene=ensembl[ids], GeneName=genename[ids], GeneSymbol=genesymbol[ids])





#parse the genomic locations. Some span exon junctions and so we'll put
#their 2 genomic locations in as separate entries in the RangedData object...

parse.genomic.ranges <- function(dat){
  g.ranges <- dat[,"Genomic_location"]
  g.ranges <- sapply(dat[,"Genomic_location"], function(x){strsplit(x,"[, ]")})
  g.ranges <- sapply(g.ranges, function(x){strsplit(x," ")})
  split.ranges.lengths <- sapply(g.ranges, function(x){length(x)})
  split.ranges.inds <- rep(1:length(g.ranges), split.ranges.lengths)
  dat <- dat[split.ranges.inds,]
  g.ranges <- unlist(g.ranges)
  g.ranges <- sapply(g.ranges, function(x){strsplit(x,":")})
  g.ranges <- do.call("rbind",g.ranges) 
  colnames(g.ranges) <- qw(Chr, Start, End, Strand)
  top.strand <- which(g.ranges[,"Strand"]=="+")
  bottom.strand <- which(g.ranges[,"Strand"]=="-")
  g.ranges[top.strand,"Strand"] <- 1
  g.ranges[bottom.strand, "Strand"] <- -1
  dat <- dat[,colnames(dat)!="Genomic_location"]
  return(list(dat,g.ranges))
}

f.test.annot.granges <- parse.genomic.ranges(f.test.annot)

t.test.1.annot.granges <- parse.genomic.ranges(t.test.1.annot)
t.test.2.annot.granges <- parse.genomic.ranges(t.test.2.annot)
t.test.3.annot.granges <- parse.genomic.ranges(t.test.3.annot)


#create a RangedData object for the region start and end:
library(IRanges)


make.rd <- function(dat.granges, filename){
  dat <- dat.granges[[1]]
  g.ranges <- dat.granges[[2]]
  
  nms <- paste(g.ranges[,"Chr"], paste(g.ranges[,"Start"],g.ranges[,"End"], sep="-"), g.ranges[,"Strand"], sep=":")

  ## there are a couple of instances where different (but very similar) probes are mapping
  ## to the exact same region of the genome. So we can't just use the position as a name:
  #wtf <- nms[which(duplicated(nms))]
  #wtf <- which(nms %in% wtf)
  #limma.annot[wtf,]
  #g.ranges[wtf]

  nms <- paste(nms," (", dat[,"IlluminaID"] ,")", sep="")

  #we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
  rd <- RangedData(ranges = IRanges(
                     start= as.numeric(g.ranges[,"Start"]),
                     end = as.numeric(g.ranges[,"End"]),
                     names = nms,
                     ),
                   space = g.ranges[,"Chr"],
                   values = cbind(dat, g.ranges[,"Strand"]),
                   universe = "mm9"
                   )

  # save the results as RangedData and csv
  save(rd, file=paste(filename, ".R", sep=""))

  rownames(g.ranges) <- rownames(dat) <- nms
  write.csv(cbind(g.ranges,dat), file=paste(filename, ".csv",sep=""), row.names=F)
}


make.rd(f.test.annot.granges, 'f_test_limma_rd')
make.rd(t.test.1.annot.granges, 't_test_bmpvctx_limma_rd')
make.rd(t.test.2.annot.granges, 't_test_fbsvctx_limma_rd')
make.rd(t.test.3.annot.granges, 't_test_bmpvfbs_limma_rd')
       

