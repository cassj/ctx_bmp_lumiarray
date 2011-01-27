#!/usr/bin/Rscript

qw <- function(...) {
  as.character(sys.call()[-1])
}

options(stringsAsFactors = FALSE);
args <- commandArgs(trailingOnly=TRUE)

# load datafile from pp_expression_data
limmafile <- args[1]
remoat.annot.file <- args[2]

#limmafile <- "/mnt/data/limma_results.csv"
#remoat.annot.file <- "/mnt/work/lib/Annotation_Illumina_Mouse-WG-V1_mm9_V1.0.0_Aug09.txt"

limma <- read.csv(limmafile)
limma<-limma[,-1]
colnames(limma)[1]<-"IlluminaID"






######
# Annotation from ReMOAT

# BTW, if we assume that the ~46K probeIDs in the remoat.annot correspond to the entire array
# and that most of those have a 1:1 relationship to a targetID then the "raw" data from KY
# has been seriously filtered. The entire limma table contains only 15060 genes.
# Nothing we can do about it unless we can get hold of the earlier data. Just be wary when
# interpreting results I guess...

# mappings to mm9

# Annotate from the ReMOAT data
# see http://nar.oxfordjournals.org/cgi/content/full/gkp942
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
limma.annot = merge(limma, remoat.annot, by.x="IlluminaID", by.y="Probe_id")

#ok, we actually have probeIds so the relationship should be 1:1

# We probably don't need *all* the annotation, although we do need to have the genome and
# transcriptome position of the probe.
# Multiple definitions of transcriptome though, so we need results for all of them.
# see http://www.compbio.group.cam.ac.uk/Resources/Annotation/Description_1.0.0.xls
# for colname definitions

nms <- qw(IlluminaID, logFC, AveExpr, t, P.Value, adj.P.Val, B,Lumi_id,Probe_sequence,Probe_type, Quality_score, Genomic_location,
          RefSeq_transcripts, Proportion_RefSeq_transcripts,
          UCSC_transcripts,Proportion_UCSC_transcripts,
          GenBank_transcripts, Proportion_GenBank_transcripts,
          Exons,
          Ensembl_transcripts, Proportion_Ensembl_transcripts,
          Lumi_transcriptomic_annotation, Lumi_transcriptomic_match )

limma.annot <- limma.annot[,nms]


#parse the genomic locations. Some span exon junctions and so we'll put
#their 2 genomic locations in as separate entries in the RangedData object...
g.ranges <- limma.annot[,"Genomic_location"]
g.ranges <- sapply(limma.annot[,"Genomic_location"], function(x){strsplit(x,",")})


# split up the ranges where necessary
split.ranges.lengths <- sapply(g.ranges, function(x){length(x)})

# 14 of these have no genomic range - although they hit a transcript.
# None of them look very interesting (low expression, no evidence for \de)
# so I'll just drop them
split.ranges.inds <- rep(1:length(g.ranges), split.ranges.lengths)

# duplicate the rows in the data table where we have split ranges
limma.annot <- limma.annot[split.ranges.inds,]
g.ranges <- unlist(g.ranges)

# parse the g.ranges into "chr","start","end","strand"
g.ranges <- sapply(g.ranges, function(x){strsplit(x,":")})
g.ranges <- do.call("rbind",g.ranges) 
colnames(g.ranges) <- qw(Chr, Start, End, Strand)
top.strand <- which(g.ranges[,"Strand"]=="+")
bottom.strand <- which(g.ranges[,"Strand"]=="-")
g.ranges[top.strand,"Strand"] <- 1
g.ranges[bottom.strand, "Strand"] <- -1

#get rid of the old genomic pos
limma.annot <- limma.annot[,colnames(limma.annot)!=("Genomic_location")]


#create a RangedData object for the region start and end:
library(IRanges)


nms <- paste(g.ranges[,"Chr"], paste(g.ranges[,"Start"],g.ranges[,"End"], sep="-"), g.ranges[,"Strand"], sep=":")

## there are a couple of instances where different (but very similar) probes are mapping
## to the exact same region of the genome. So we can't just use the position as a name:
#wtf <- nms[which(duplicated(nms))]
#wtf <- which(nms %in% wtf)
#limma.annot[wtf,]
#g.ranges[wtf]

nms <- paste(nms," (", limma.annot[,"IlluminaID"] ,")", sep="")

#we have to give them names to avoid a bug in ChIPpeakAnnot if we want to use it later
rd.limma <- RangedData(ranges = IRanges(
                         start= as.numeric(g.ranges[,"Start"]),
                         end = as.numeric(g.ranges[,"End"]),
                         names = nms,
                         ),
                       space = g.ranges[,"Chr"],
                       values = cbind(limma.annot, g.ranges[,"Strand"]),
                       universe = "mm9"
                       )

# save the results as RangedData and csv
save(rd.limma, file="limma_rd.R")

rownames(g.ranges) <- rownames(limma.annot) <- nms
write.csv(cbind(g.ranges,limma.annot), file="limma_rd.csv", row.names=F)



