
typemlst <- function(genomePath,type)
{

  if(type=='ecoli-warwick'){genesPath <- list.files(system.file("extdata/typing/Escherichia-coli/MSLT-warwick/",package = "Pathogenomics"),full.names=T) }
  if(type=='klebsiella-pasteur'){genesPath <- list.files(system.file("extdata/typing/Klebsiella-pneumoniae/MSLT-pasteur/",package = "Pathogenomics"),full.names=T) }

  combination <- genesPath[grep('txt',genesPath)]

  genesPath <- genesPath[grep('fasta',genesPath)]
  genesName  <- gsub(pattern='.fasta',replacement='',x=basename(genesPath))

  try(unlink("temp", recursive=TRUE))
  dir.create('temp')
  dir.create('temp/dbblast')
  myarg <-paste0('-in ',genomePath,' -out temp/dbblast/db -dbtype nucl')
  system2(command = 'makeblastdb', args = myarg,stdout=F)

  Ngene <- length(genesPath)
  result <- character()

  for (i in 1:Ngene)
  {
    myarg <-  paste0('-query ',genesPath[i],' -db temp/dbblast/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc qlen length pident qstart qend sacc sstart send "' )
    system2(command = 'blastn', args = myarg)

    blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
    if(class(blast)=='data.frame')
    {
      colnames(blast) <- c('querry.access','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
      blast <- blast[(blast$alignment.lenght==blast$querry.length)&(blast$pc.ident.==100),]
      result[i] <- paste0('',as.character(blast$querry.access))
    }
    else{result[i] <- ''}
  }

  combination <- read.table(combination,sep='\t',header=T)

  result
  resultnum <- substr(result,start=unlist(lapply(gregexpr("[0123456789]",result), function(l) l[[1]])),stop=nchar(result))
  resultchar <- substr(result,start=1,stop=unlist(lapply(gregexpr("[0123456789]",result), function(l) l[[1]]))-1)

  idx <- numeric()
  for(i in 1:Ngene)
    {
    idx[i] <- grep(resultchar[i],colnames(combination))
    }

  combinationReduced <- combination[,idx]
  rownames(combinationReduced) <- combination$ST

  result <- combinationReduced[apply(mapply("==",as.numeric(resultnum),combinationReduced),1,sum)==Ngene,]
  result$ST <- rownames(result)
  rownames(result) <- NULL
  return(result)

}



