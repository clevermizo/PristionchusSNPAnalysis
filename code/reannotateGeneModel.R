reannotateGeneModel = function(gdata,k=2000){
gndata=gdata
ntranscript = length(unique(gndata$transcript_id))
transcriptnames = unique(gndata$transcript_id)
tmp = data.frame(seqid = rep(gndata$seqid[1],ntranscript),
                 source = rep(gndata$source[1],ntranscript),
                 feature = rep("promoter",ntranscript),
                 start = rep(NA,ntranscript),
                 end = rep(NA,ntranscript),
                 score = rep(gndata$score[1],ntranscript),
                 strand = rep(gndata$strand[1],ntranscript),
                 frame = rep(gndata$frame[1],ntranscript),
                 transcript_id = transcriptnames,
                 gene_id = rep(gndata$gene_id[1],ntranscript),
                 gene_name = rep(gndata$gene_name[1],ntranscript))
for(j in 1:ntranscript){
  tmp2 = gndata[gndata$transcript_id == unique(gndata$transcript_id)[j],]
  if( tmp2$strand[1] == "+" ){
    tmp$start[j] = min(tmp2$start) - k
    tmp$end[j] = min(tmp2$start) - 1
  }else if( tmp2$strand[1] == "-"){
    tmp$end[j] = max(tmp2$end) + k
    tmp$start[j] = max(tmp2$end) + 1
  }
}
gndata = rbind(gndata,tmp)
gndata = gndata[order(gndata$transcript_id,gndata$seqid,gndata$start),]
cdsdata = gndata[gndata$feature=="CDS",]
gndata = gndata[gndata$feature != "CDS",]
gndata_transcripts = vector("list",length = length(unique(gndata$transcript_id)))
names(gndata_transcripts) = transcriptnames

for(j in transcriptnames){
  tmp = gndata[gndata$transcript_id==j,]
  if(tmp$start[tmp$feature=="promoter"]<0){
    tmp$start[tmp$feature=="promoter"]=1
  }
  cds = cdsdata[cdsdata$transcript_id==j,]
  
  exons = tmp[tmp$feature=="exon",]
  if(nrow(exons)>1){
  introns = exons
  introns$start = NA
  introns$end = NA
  introns$feature = "intron"
  introns = introns[1:(nrow(introns)-1),]
  introns$start = exons$end[1:(nrow(exons)-1)]+1
  introns$end = exons$start[2:nrow(exons)]-1
  }
  
  if(tmp$strand[1] == "+"){
    cdss = cds[1,]
    cdsstop = cds[nrow(cds),]
    
    tmp2 = tmp
    cdss_test = cdss$start - tmp$start
    if(sum(cdss_test==0)==0){
      cdss_bin = which(cdss_test<0)[1]-1
      tmp2$end[cdss_bin] = cdss$start - 1
      tmp2$feature[2:cdss_bin] = "5utr"
      tmp2 = rbind(tmp2,cdss)
      tmp2 = tmp2[order(tmp2$seqid,tmp2$start),]
    }else{
      cdss_bin = which(cdss_test==0)
      tmp2$feature[cdss_bin]="CDS"
      if("exon" %in% tmp2$feature[1:cdss_bin]){
      tmp2$feature[2:(cdss_bin-1)]="5utr"
      }
    }
    
    cdsstop_test = cdsstop$end - tmp2$end
    if(sum(cdsstop_test==0)==0){
      cdsstop_bin = which(cdsstop_test < 0)[1]
      tmp2$start[cdsstop_bin] = cdsstop$end + 1
      tmp2$feature[cdsstop_bin:nrow(tmp2)] = "3utr"
      tmp2 = rbind(tmp2,cdsstop)
      tmp2 = tmp2[order(tmp2$seqid,tmp2$start),]
    }else{
      cdsstop_bin = which(cdsstop_test ==0)
      tmp2$feature[cdsstop_bin] = "CDS"
      if("exon" %in% tmp2$feature[cdsstop_bin:nrow(tmp2)]){
        tmp2$feature[(cdsstop_bin+1):nrow(tmp2)]="3utr"
      }
    }
    cdsidx = which(tmp2$feature=="CDS")
    if(length(cdsidx)>1){
    tmp2$feature[cdsidx[1]:cdsidx[2]]="CDS"
    }
    
  }else if(tmp$strand[1] == "-"){
    cdss = cds[nrow(cds),]
    cdsstop = cds[1,]
    tmp2=tmp
    cdss_test = cdss$end - tmp$end
    if(sum(cdss_test==0)==0){
      cdss_bin = which(cdss_test<0)[1]
      tmp2$start[cdss_bin] = cdss$end+1
      tmp2$feature[cdss_bin:(nrow(tmp2)-1)] = "5utr"
      tmp2 = rbind(tmp2,cdss)
      tmp2 = tmp2[order(tmp2$seqid,tmp2$start),]
    }else{
      cdss_bin = which(cdss_test==0)
      tmp2$feature[cdss_bin]="CDS"
      if("exon" %in% tmp2$feature[cdss_bin:nrow(tmp2)]){
        tmp2$feature[(cdss_bin+1):(nrow(tmp2)-1)]="5utr"
      }
    }
    
    cdsstop_test = cdsstop$start - tmp2$start
    if(sum(cdsstop_test==0)==0){
      cdsstop_bin = which(cdsstop_test < 0)[1]-1
      tmp2$end[cdsstop_bin] = cdsstop$start - 1
      tmp2$feature[1:cdsstop_bin] = "3utr"
      tmp2 = rbind(tmp2,cdsstop)
      tmp2 = tmp2[order(tmp2$seqid,tmp2$start),]
    }else{
      cdsstop_bin = which(cdsstop_test ==0)
      tmp2$feature[cdsstop_bin] = "CDS"
      if("exon" %in% tmp2$feature[1:cdsstop_bin]){
        tmp2$feature[1:(cdsstop_bin-1)]="3utr"
      }
    }
    cdsidx = which(tmp2$feature=="CDS")
    if(length(cdsidx)>1){
    tmp2$feature[cdsidx[1]:cdsidx[2]]="CDS"
    }
  }
  
  if(nrow(exons)>1){
  tmp2 = rbind(tmp2,introns)
  }
  tmp2 = tmp2[order(tmp2$seqid,tmp2$start),]
  gndata_transcripts[[j]] = tmp2
}

return(gndata_transcripts)
}



