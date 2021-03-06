extractgenbank <- function(minlength,maxlength,directory,myquery)
{
  library('seqinr')
  library('ape')
  library('Biostrings')

  choosebank("genbank")
  fc <- query("fc",myquery)

  names <- unlist(fc$req)
  Nseq <- fc$nelem
  mylength <- numeric()
  for(i in 1:Nseq){mylength[i] <- attr(fc$req[[i]],"length")}

  tokeep <- (mylength>minlength)&(mylength<maxlength)

  indextokeep <- which(tokeep)
  Nseq <- length(indextokeep)
  sequences <- DNAStringSet()

  print(Nseq)
  for (j in 1:Nseq)
  {
    sequences <- c(sequences, DNAStringSet(getSequence(fc$req[[indextokeep[j]]], as.string = TRUE)))
    print(j)
  }

  annotation <- list()
  for (j in 1:Nseq)
  {
    annotation1seq <- getAnnot(fc$req[[indextokeep[j]]],nbl = 15)
    annotation1seq <- annotation1seq[grep('DEFINITION',annotation1seq)]
    annotation[[j]] <- annotation1seq
    print(j)
  }

  annotation <- unlist(annotation)
  annotation <- gsub(pattern=', complete genome.',replacement='',annotation)
  annotation <- gsub(pattern=', complete sequence.',replacement='',annotation)
  annotation <- gsub(pattern='DEFINITION  ',replacement='',annotation)
  annotation <- gsub(pattern='/',replacement='',annotation)


  names <- names[indextokeep]
  names <- paste(names,annotation)
  names <- gsub(pattern=' ', replacement ='_' ,names)
  names <- gsub("[[:punct:]]", ".", names)

  names(sequences) <- names

  for(i in 1:Nseq)
  {
    writeXStringSet(sequences[i],paste0(directory,names[i],'.fasta'),format='fasta')
  }
}
