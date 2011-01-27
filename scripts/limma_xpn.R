#!/usr/bin/Rscript

library(beadarray)

# because I can't be bothered to type all the quotes
qw <- function(...) {
  as.character(sys.call()[-1])
}

options(stringsAsFactors = FALSE);
options(scipen=10)

args <- commandArgs(trailingOnly=TRUE)
filename = args[1]
outfile = args[2]


filename<-"/mnt/data/Non-background_subtracted_Sample_Probe_Profile.txt"
outfile <- "/mnt/data/limma_results.csv"

test <- read.csv(filename, skip=7)

BSData<-readBeadSummaryData(filename,
                            sep="\t",
                            skip=7,
                            ProbeID="ProbeID",
                            columns=list(
                              exprs="AVG_Signal",
                              se.exprs="BEAD_STDEV",
                              NoBeads="Avg_NBEADS",
                              Detection="Detection"
		   )
)	

E = normaliseIllumina(BSData, method="quantile", transform="log2")
data <- exprs(E)

library(limma)

chips <- colnames(exprs(BSData))

BMP4 <- c("4203884026_A", "4203884026_B", "4203884035_A", "4203884035_B")
FBS <- c("4203884026_C", "4203884026_D", "4203884035_C", "4203884035_D")
CTX12 <- c("4203884026_E", "4203884026_F", "4203884035_E", "4203884035_F")

bmp <- c(1,1,rep(0,4),1,1,rep(0,4))
fbs <- c(0,0,1,1,rep(0,4),1,1,0,0)
ctx <- c(rep(0,4),1,1,rep(0,4),1,1)

design<-cbind(bmp, fbs, ctx)
rownames(design) <- chips


fit<-lmFit(exprs(E), design)
cont.matrix<-makeContrasts(bmpvctx=bmp-ctx,
                           fbsvctx=fbs-ctx,
                           bmpvfbs=bmp-fbs,
                           levels=design)
fit<-contrasts.fit(fit, cont.matrix)
ebFit<-eBayes(fit)


#F test
#for now, consider the t-tests on their own when correcting. Unsure if this is best or not.

#Wait - we're not actually using this?
#write.fit(ebFit, file="limma_fit.csv", F.adjust="BH", adjust="BH", method="separate", sep=",")
#data<-read.table(outfile, sep="\t", header=T)

data<- topTable(ebFit, number=nrow(data))
write.csv(data,outfile)

#Individual comparisons:
outfile.1 <- sub('.csv', '_bmp_v_ctx.csv', outfile) 
outfile.2 <- sub('.csv', '_fbs_v_ctx.csv', outfile) 
outfile.3 <- sub('.csv', '_bmp_v_fbs.csv', outfile) 

data<- topTable(ebFit, coef=1, number=nrow(data))
write.csv(data,outfile.1)

data<- topTable(ebFit, coef=2, number=nrow(data))
write.csv(data,outfile.2)

data<- topTable(ebFit, coef=3, number=nrow(data))
write.csv(data,outfile.3)










