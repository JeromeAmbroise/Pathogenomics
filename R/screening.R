#' Screen for a set of genes of interst (e.g. virulence genes) from a (draft of complete) genome
#'
#' screening returns a vector indicating if the genome matches each gene of interst.
#'
#' @param genomePath A string of character corresponding to the path of the FASTA file of the genome to be screened.
#'
#' @param genesPath A vector of string of character corresponding to the paths of the FASTA files of the genes of interst (e.g. virulence genes).
#'
#' @param lengthconf The minimum accepted length of the hit (in percentage of the gene of interest).
#'
#' @param identconf The minimum accepted identity of the hit (in percentage of identity).
#' @param outputdir The directory where output files are stored.
#'
#' @return result is a vector containing 0 and 1. This vector has the length of the number of genes of interest.
#' The name of the elements of this vector corresponds to the name of the gene of interest.
#' @examples
#' genomePath <- list.files(system.file("extdata/genome", package = "Pathogenomics"),full.names=T)
#' Ngenomes <- length(genomePath)
#' genesPath <- list.files(system.file("extdata/virulence/Escherichia-coli/Virulencefinder", package = "Pathogenomics"),full.names=T)
#' Ngenes <- length(genesPath)
#' screening(genomePath= genomePath[1],genesPath=genesPath,lengthconf = 60,identconf = 90,outputdir = '7-virulence-Irenge')

#' @export
screening <- function(genomePath,genesPath,lengthconf,identconf,outputdir)
{
#Add comments and doc
  genomeName <- gsub(pattern='.fasta',replacement='',x=basename(genomePath))
  genesName  <- gsub(pattern='.fasta',replacement='',x=basename(genesPath))

  if(dir.exists(outputdir)==F){dir.create(outputdir)}

  try(unlink("temp", recursive=TRUE))
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  Ngene <- length(genesPath)
  result <- numeric()

  for (i in 1:Ngene)
  {
    myarg <-  paste0('-query ',genesPath[i],' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc bitscore qlen length pident qstart qend sacc sstart send "' )
    system2(command = 'blastn', args = myarg)

    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    if (class(blast) == 'data.frame')
    {
      colnames(blast) <- c('querry.access','bitscore','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      blast <- blast[sort.list(blast$bitscore, decreasing = T), ]
      pc.length <- 100 * round(blast$alignment.lenght / blast$querry.length, 3)
      blast <-  data.frame(blast[, c(1, 2, 3)], pc.length, blast[, -c(1, 2, 3)])
      blast <- blast[blast$pc.length >= lengthconf, ]
      blast <- blast[blast$pc.ident. >= identconf, ]
      if (dim(blast)[1] > 0)
      {
        #write.csv(blast,paste0(outputdir,'/',paste0(genomeName,genesName[i]),'.csv'),row.names = F)
        result[i] <- blast$pc.length[1]
      }
      else(result[i] <- 0)
    }
    else(result[i] <- 0)
    try(file.remove("temp/blast.txt"))
  }
  names(result) <- genesName

  try(unlink("temp", recursive=TRUE))
  return(result)
}









findlocation <- function(genomePath,genesPath,lengthconf,identconf,outputdir)
{
  #Add comments and doc
  genomeName <- gsub(pattern='.fasta',replacement='',x=basename(genomePath))
  genesName  <- gsub(pattern='.fasta',replacement='',x=basename(genesPath))

  if(dir.exists(outputdir)==F){dir.create(outputdir)}

  try(unlink("temp", recursive=TRUE))
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  Ngene <- length(genesPath)
  result <- character()

  for (i in 1:Ngene)
  {
    myarg <-  paste0('-query ',genesPath[i],' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc bitscore qlen length pident qstart qend sacc sstart send "' )
    system2(command = 'blastn', args = myarg)

    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    if (class(blast) == 'data.frame')
    {
      colnames(blast) <- c('querry.access','bitscore','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      blast <- blast[sort.list(blast$bitscore, decreasing = T), ]
      pc.length <- 100 * round(blast$alignment.lenght / blast$querry.length, 3)
      blast <-  data.frame(blast[, c(1, 2, 3)], pc.length, blast[, -c(1, 2, 3)])
      blast <- blast[blast$pc.length >= lengthconf, ]
      blast <- blast[blast$pc.ident. >= identconf, ]
      blast <- blast[sort.list(blast$bitscore,decreasing=T), ]

      if (dim(blast)[1] > 0)
      {
        write.csv(blast,paste0(outputdir,'/',paste0(genomeName,genesName[i]),'.csv'),row.names = F)
        result[i] <- paste0('% ident = ',blast$pc.ident.[1],'; % cover. = ',blast$pc.length[1],'; contig = ',blast$subject.access[1],'; start = ',blast$subject.start[1],'; end = ',blast$subject.end[1])
      }
      else(result[i] <- 0)
    }
    else(result[i] <- 0)
    try(file.remove("temp/blast.txt"))
  }
  names(result) <- genesName

  try(unlink("temp", recursive=TRUE))
  return(result)
}










find.closest <- function(genomePath,genepath,lengthconf = 100, identconf =100 )
{
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  myarg <-  paste0('-query ',genepath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc bitscore qlen length pident qstart qend sacc sstart send "' )
  system2(command = 'blastn', args = myarg)

  blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
  if(class(blast)=='data.frame')
  {
    colnames(blast) <- c('querry.access','bit.score','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
    blast <- blast[(blast$alignment.lenght>=blast$querry.length*lengthconf/100),]
    blast <- blast[(blast$pc.ident.>=identconf),]
    result <- paste0('',as.character(blast$querry.access)[which.max(blast$bit.score)])
  }
  else{result <- ''}
  return(result)

}





extract.closest <- function(genomePath,genepath,lengthconf = 95, identconf =95,offset=0 )
{
  library(Biostrings)
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  myarg <-  paste0('-query ',genepath,' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "' )
  system2(command = 'blastn', args = myarg)

  blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
  if(class(blast)=='data.frame')
  {
    colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
    blast <- blast[(blast$alignment.lenght>=lengthconf/100*blast$querry.length),]
    blast <- blast[(blast$pc.ident.>=identconf),]
    blast <- blast[sort.list(blast$alignment.lenght,decreasing=T),]
    if(dim(blast)[1]>=1)
    {
      sequence <- readDNAStringSet(genomePath)
      sequence <- sequence[names(sequence)==as.character(blast$subject.access[1])]


      if(blast$subject.start[1]<blast$subject.end[1])
      {sequence <- subseq(sequence,start=max(1,blast$subject.start[1]-offset),end=min(blast$subject.end[1]+offset,width(sequence)))}

      if(blast$subject.start[1]>blast$subject.end[1])
      {sequence <- subseq(sequence,start=max(1,blast$subject.end[1]-offset),end=min(blast$subject.start[1]+offset,width(sequence)))
      sequence <- reverseComplement(sequence)}

    }
    else{sequence <- ''}
  }
  else{sequence <- ''}
  sequence <- DNAStringSet(sequence)
  return(sequence)
}










screeningp <- function(genomePath,genesPath,lengthconf,identconf,outputdir)
{
  #Add comments and doc
  genomeName <- gsub(pattern='.fasta',replacement='',x=basename(genomePath))
  genesName  <- gsub(pattern='.fasta',replacement='',x=basename(genesPath))

  if(dir.exists(outputdir)==F){dir.create(outputdir)}

  try(unlink("temp", recursive=TRUE))
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype prot')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  Ngene <- length(genesPath)
  result <- numeric()

  for (i in 1:Ngene)
  {
    myarg <-  paste0('-query ',genesPath[i],' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc bitscore qlen length pident qstart qend sacc sstart send "' )
    system2(command = 'blastp', args = myarg)

    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    if (class(blast) == 'data.frame')
    {
      colnames(blast) <- c('querry.access','bitscore','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      blast <- blast[sort.list(blast$bitscore, decreasing = T), ]
      pc.length <- 100 * round(blast$alignment.lenght / blast$querry.length, 3)
      blast <-  data.frame(blast[, c(1, 2, 3)], pc.length, blast[, -c(1, 2, 3)])
      blast <- blast[blast$pc.length >= lengthconf, ]
      blast <- blast[blast$pc.ident. >= identconf, ]
      if (dim(blast)[1] > 0)
      {
        write.csv(blast,paste0(outputdir,'/',paste0(genomeName,genesName[i]),'.csv'),row.names = F)
        result[i] <- 1
      }
      else(result[i] <- 0)
    }
    else(result[i] <- 0)
    try(file.remove("temp/blast.txt"))
  }
  names(result) <- genesName

  try(unlink("temp", recursive=TRUE))
  return(result)
}


