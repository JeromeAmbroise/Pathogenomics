typemlst <- function(genomePath,type)
{

  if(type=='ecoli-warwick'){genesPath <- list.files(system.file("extdata/typing/Escherichia-coli/MSLT-warwick/",package = "Pathogenomics"),full.names=T) }
  if(type=='klebsiella-pasteur'){genesPath <- list.files(system.file("extdata/typing/Klebsiella-pneumoniae/MLST-pasteur/",package = "Pathogenomics"),full.names=T) }
  if(type=='vibrio-pasteur'){genesPath <- list.files(system.file("extdata/typing/Vibrio-cholerae/MLST-pasteur/",package = "Pathogenomics"),full.names=T) }
  if(type=='bacillus-th'){genesPath <- list.files(system.file("extdata/typing/Bacillus/MLST-Tourasse-Helgason",package = "Pathogenomics"),full.names=T) }
  if(type=='bacillus-pasteur'){genesPath <- list.files(system.file("extdata/typing/Bacillus/MLST-pasteur",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncACcgPMLST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncACcgPMLST/",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncACPMLST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncACPMLST/",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncFRST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncFRST/",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncHI1MLST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncHI1MLST/",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncHI2DLST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncHI2DLST/",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncI1MLST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncI1MLST/",package = "Pathogenomics"),full.names=T) }
  if(type=='plasmid-IncNMLST'){genesPath <- list.files(system.file("extdata/typing/plasmid/IncNMLST/",package = "Pathogenomics"),full.names=T) }


  if(exists('genesPath')==F)
  {cat('specify one of the following type: ecoli-warwick vibrio-pasteur klebsiella-pasteur plasmid-IncACcgPMLST plasmid-IncACPMLST plasmid-IncFRST plasmid-IncHI1MLST plasmid-IncHI2DLST plasmid-IncI1MLST plasmid-IncNMLST')  }
  else
  {
    combination <- genesPath[grep('txt',genesPath)]

    genesPath <- genesPath[grep('fas',genesPath)]
    genesName  <- gsub(pattern='.fasta',replacement='',x=basename(genesPath))
    genesName  <- gsub(pattern='.fas',replacement='',x=basename(genesName))

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
        result[i] <- paste0('',as.character(blast$querry.access[which.max(blast$querry.length)]))
      }
      else{result[i] <- ''}
    }

    combination <- read.table(combination,sep='\t',header=T)
    combination <- combination[duplicated(combination$ST)==F,]
    resultnum <- substr(result,start=unlist(lapply(gregexpr("[0123456789]",result), function(l) l[[1]])),stop=nchar(result))
    names(resultnum) <- genesName

    if(min(nchar(resultnum))>0)
    {
      idx <- numeric()
      for(i in 1:Ngene){idx[i] <- grep(genesName[i],colnames(combination))}
      combinationReduced <- combination[,idx]
      rownames(combinationReduced) <- combination$ST
      ST <- paste0('',rownames(combinationReduced[apply(mapply("==",as.numeric(resultnum),combinationReduced),1,sum)==Ngene,]))
    }
    else{ST<-''}
    names(ST) <- 'ST'

    result <- c(resultnum,ST)
    return(result)

  }
}
