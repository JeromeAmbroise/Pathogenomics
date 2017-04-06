typemlst <- function(genomePath,type)
{

  #### ESCHERICHIA COLI

  if(type=='ecoli-warwick')
  {
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/ADK.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    ADK <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      ADK <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==536)&(blast$pc.ident.==100)]))[1]
      print('ADK done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/FUMC.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    FUMC <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      FUMC <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==469)&(blast$pc.ident.==100)]))[1]
      print('FUMC done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/GYRB.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    GYRB <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      GYRB <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==460)&(blast$pc.ident.==100)]))[1]
      print('GYRB done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/ICD.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    ICD <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      ICD <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==518)&(blast$pc.ident.==100)]))[1]
      print('ICD done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/MDH.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    MDH <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      MDH <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==452)&(blast$pc.ident.==100)]))[1]
      print('MDH done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/PURA.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    PURA <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      PURA <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==478)&(blast$pc.ident.==100)]))[1]
      print('PURA done')
    }



    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Escherichia-coli/MSLT-warwick/RECA.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    RECA <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      RECA <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==510)&(blast$pc.ident.==100)]))[1]
      print('RECA done')
    }

    try(unlink("temp", recursive=TRUE),silent=T)
    MSLT <- c(ADK,FUMC,GYRB,ICD,MDH,PURA,RECA)
    return(MSLT)
  }


  #### KLEBSIELLA PASTEUR

  if(type=='klebsiella-pasteur')
  {

    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/GAPA.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    GAPA <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      GAPA <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==450)&(blast$pc.ident.==100)]))[1]
      print('GAPA done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/INFB.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    INFB <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      INFB <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==318)&(blast$pc.ident.==100)]))[1]
      print('INFB done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/MDH.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    MDH <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      MDH <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==477)&(blast$pc.ident.==100)]))[1]
      print('MDH done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/PGI.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    PGI <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      PGI <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==432)&(blast$pc.ident.==100)]))[1]
      print('PGI done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/PHOE.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    PHOE <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      PHOE <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==420)&(blast$pc.ident.==100)]))[1]
      print('PHOE done')
    }


    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/RPOB.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    RPOB <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      RPOB <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==501)&(blast$pc.ident.==100)]))[1]
      print('RPOB done')
    }



    try(rm(blast),silent=T)
    try(unlink("temp", recursive=TRUE),silent=T)
    myarg <-paste0('-in ',system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/TONB.fasta", package = "Pathogenomics"),' -out temp/dbblast/db -dbtype nucl')
    system2(command = 'makeblastdb', args = myarg,stdout=F)
    myarg <- paste0('-query ',genomePath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "')
    system2(command='blastn',args=myarg)
    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    TONB <- ''
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      TONB <- paste0('',as.character(blast$subject.access[(blast$alignment.lenght==414)&(blast$pc.ident.==100)]))[1]
      print('TONB done')
    }

    try(unlink("temp", recursive=TRUE),silent=T)
    MSLT <- c(GAPA,INFB,MDH,PGI,PHOE,RPOB,TONB)
    return(MSLT)
  }






}
