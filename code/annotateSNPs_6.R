library(stringr)
library(Biostrings)
library(Rsamtools)
options(stringsAsFactors = FALSE)
#setwd("~/mnlscdata/mrieger/RNAseq/PristionchusSNPAnalysis/code")
#### files ####
vcf="../Ppa_alignments_variants_qgeq20.vcf"
gtf="../../../Genomes/Annotation/pristionchus_pacificus.PRJNA12644.WBPS11.annotations.gtf"
genome="../../../Genomes/Sources/pristionchus_pacificus.PRJNA12644.WBPS11.genomic.fa"
k=2000

#### setup GTF ####
#gtfx = read.table(gtf,header = FALSE,sep="\t")
#gtfnames=c("seqid","source","feature","start","end","score","strand","frame","transcript_id","gene_id","gene_name")
#v9=gtfx[,9]
#gtfx = gtfx[,1:8]
#v9 = str_split(v9,";",simplify=TRUE)
#gtfx = cbind(gtfx,str_split(v9[,1]," ",simplify = TRUE)[,2])
#gtfx = cbind(gtfx,str_split(v9[,2]," ",simplify = TRUE)[,3])
#gtfx = cbind(gtfx,str_split(v9[,3]," ",simplify = TRUE)[,3])
#colnames(gtfx) = gtfnames
#rm(v9)
#saveRDS(gtfx, "pristionchus_pacificus.PRJNA12644.WBPS11.annotations.gtf.rds")
gtfx=readRDS("pristionchus_pacificus.PRJNA12644.WBPS11.annotations.gtf.rds")


#### Annotation Round 1 ####
# We are going to annotate vcf
# Categories: promoter, 3utr, 5utr, nongenic, CDS, intron

geneNames = unique(gtfx$gene_name)
#vcfdata = read.table("../Ppa_alignments_variants_qgeq20.vcf",header = TRUE,comment.char = '#',sep="\t")
#vcfdata$snpclass = "nongenic"
#vcfdata$gene_id = "nongenic"

##### Relabel Samples ######
#colnames(vcfdata)[c(10,11)] = c("rs5194_1","rs5194_2")
#vcfdata$genotype_1 = str_split(vcfdata$rs5194_1,":",simplify=TRUE)[,1]
#vcfdata$genotype_1[vcfdata$genotype_1 == "./."] = NA
#vcfdata$genotype_2 = str_split(vcfdata$rs5194_2,":",simplify=TRUE)[,1]
#vcfdata$genotype_2[vcfdata$genotype_2 == "./."] = NA
# 847,314 SNPs
# Requiring a call, with Quality >=20, for BOTH samples:
#vcfdata = vcfdata[!is.na(vcfdata$genotype_1) & !is.na(vcfdata$genotype_2),]
# 528,468 remaining SNPs

#saveRDS(vcfdata,"vcfdata_unannotated.rds")
##### load function for generating gene/transcript models #####
source("reannotateGeneModel.R")

#### Loop over Genes ####

#for(g in geneNames){
#  print(paste0("Processing SNPs in gene ",g))
#  gdata = gtfx[gtfx$gene_id==g,]
#  gndata_transcripts = reannotateGeneModel(gdata,k=2000)
#  gndata_collapse = do.call(rbind,gndata_transcripts)
#  idx =  which(vcfdata$CHROM == gdata$seqid[1] & vcfdata$POS >= min(gndata_collapse$start) &
#                   vcfdata$POS <= max(gndata_collapse$end))
#  if(length(idx)>0){
#  snps = vcfdata[idx,]
#  for(s in 1:nrow(snps)){
#    lwr = which(gndata_collapse$start <= snps$POS[s] & gndata_collapse$seqid == snps$CHROM[s])
#    tmp = gndata_collapse[lwr,]
#    uppr = which(tmp$end >= snps$POS[s] & tmp$seqid == snps$CHROM[s])
#    tmp = tmp[uppr,]
#    snps$snpclass[s] = str_c(unique(tmp$feature),collapse = ",")
#  }
#  vcfdata$snpclass[idx] = snps$snpclass
#  vcfdata$gene_id[idx] = g
#  }
#}
#saveRDS(vcfdata,"vcfdata_annotated_round1.rds")
#
vcfdata=readRDS("vcfdata_annotated_round1.rds")
#### Annotate CDS for mutation type ####
# Nonsense
# Missense
# Silent
# This requires generating a gene model, fetching the CDS by transcript
# Then, generating a second CDS with the allele altered
# Translating
# Comparing the translations.
# Save the annotation as CDS;AA1>AA2;silent|nonsense|missense
# This can be parsed and counted later, or even just grep'd.

cdsdata = vcfdata[vcfdata$snpclass == "CDS",]
cdsdata = cdsdata[50001:60000,]
rm(vcfdata)
# 77,455 SNPs in this class
for(snp in 1:nrow(cdsdata)){
    #print(paste0(snp,"/",nrow(cdsdata)))
  
    gm = reannotateGeneModel(gtfx[gtfx$gene_id == cdsdata$gene_id[snp],])
    gm.names=names(gm)
    gm.seq = vector("list",length(gm))
    names(gm.seq)=gm.names
    gm.trans = vector("list",length(gm))
    names(gm.trans)=gm.names
    gm.altid = rep(NA,length(gm))
    names(gm.altid) = gm.names
    
    matchSummary=""
    snpVal = as.character(str_split(cdsdata[snp,"ALT"],",",simplify=TRUE))
  for(alts in 1:length(snpVal)){
   for(xscript in gm.names){
    xscript.cds = gm[[xscript]][gm[[xscript]]$feature == "CDS",]
    xscript.idx = xscript.cds$start[1]:xscript.cds$end[1]
    if(nrow(xscript.cds)>1){
     for(j in 2:nrow(xscript.cds)){
        xscript.idx = c(xscript.idx,xscript.cds$start[j]:xscript.cds$end[j])
     }
    }
    if(length(xscript.idx)%%3 == 0){
        #xscript.seq = scanFa(genome,
    #                     GRanges(
    #                       seqnames = xscript.cds$seqid,
    #                       ranges = IRanges(start=xscript.cds$start,end=xscript.cds$end),
    #                       strand = xscript.cds$strand))
    # For some reason, "strand" does nothing?
    xscript.seq = scanFa(genome,
                         GRanges(
                           seqnames = xscript.cds$seqid,
                           ranges = IRanges(start=xscript.cds$start,end=xscript.cds$end)))
    xscript.seqc = str_c(as.character(xscript.seq),collapse="")
        
    if(xscript.cds$strand[1] == "-"){
     xscript.seqc = as.character(reverseComplement(DNAString(xscript.seqc)))
     xscript.idx = rev(xscript.idx)
     }
    gm.altid[xscript] = which(xscript.idx == cdsdata[snp,"POS"])    
    xscript.alt = paste0(
      substr(xscript.seqc,start=1,stop=( gm.altid[xscript]-1)),
      snpVal[alts],
      substr(xscript.seqc,start=(gm.altid[xscript]+1),stop=str_length(xscript.seqc))
      )
    
    
    gm.seq[[xscript]] = DNAStringSet(c(ref=xscript.seqc,alt=xscript.alt))
    
    gm.trans[[xscript]] = AAStringSet(
      c(ref = as.character(translate(DNAString(xscript.seqc))),
        alt = as.character(translate(DNAString(xscript.alt)))
       )
      )
    
    # Compare the two sequences
    seqmatch = mismatchTable(pairwiseAlignment(gm.trans[[xscript]]["ref"],gm.trans[[xscript]]["alt"]))
        
    # 3 cases
    
    if( nrow(seqmatch) == 0 ){
      matchSummary = c(matchSummary,paste0("CDS;silent"))
    }else{
      seqmatch = seqmatch[1,]
      if( seqmatch$SubjectSubstring == "*" ){
        matchSummary = c(matchSummary,paste0("CDS;",seqmatch$PatternSubstring,seqmatch$PatternEnd,">",seqmatch$SubjectSubstring,";nonsense"))
        }else{
        matchSummary = c(matchSummary,paste0("CDS;",seqmatch$PatternSubstring,seqmatch$PatternEnd,">",seqmatch$SubjectSubstring,";missense"))
        }
    }
    }
    }
   }
  matchSummary = matchSummary[matchSummary != ""]
  cdsdata$snpclass[snp] = str_c(unique(matchSummary),collapse="|")
  if(snp == 1){
  write.table(cdsdata[snp,],file="cdsdata_annotated.tab.6",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }else{
    write.table(cdsdata[snp,],file="cdsdata_annotated.tab.6",append=TRUE,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    }
}
#saveRDS(cdsdata,"cdsdata_annotated.rds")


