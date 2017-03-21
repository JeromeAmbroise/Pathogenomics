

Rmakeblastdb <- function(fastapath='')
{
dir.create('temp')
dir.create('temp/dbgenomes')
myarg <-paste0('-in ',fastapath,' -out temp/dbgenomes/db -dbtype nucl')
system2(command = 'makeblastdb', args = myarg)
print('DB was created')
}




screenblastdb <- function(gene,lowconflength,lowconfident,highconflength,highconfident,outputdir)
{

  myarg <-  paste0('-query ',gene,' -db temp/dbgenomes/db -out temp/blast.txt -num_threads 8 -num_alignments 10 -outfmt "7 qacc bitscore qlen length pident qstart qend sacc sstart send "' )
  system2(command = 'blastn', args = myarg)

  blast <- try(read.table('temp/blast.txt', comment.char = '#'),silent=T)
  if (class(blast) == 'data.frame')
  {
    colnames(blast) <- c('querry.access','bitscore','querry.length','alignment.lenght','pc.ident.','querry.start','querry.end','subject.access','subject.start','subject.end')
    blast <- blast[sort.list(blast$bitscore, decreasing = T), ]
    pc.length <- 100 * round(blast$alignment.lenght / blast$querry.length, 3)
    blast <-  data.frame(blast[, c(1, 2, 3)], pc.length, blast[, -c(1, 2, 3)])
    blast <- blast[blast$pc.length > lowconflength, ]
    blast <- blast[blast$pc.ident. > lowconfident, ]
    if (dim(blast)[1] > 0)
    {
      write.csv(blast,paste0(outputdir,'/',genome.name[i],'-',gene.name[j],'.csv'),row.names = F)
      if ((blast$pc.length[1]>highconflength)&(blast$pc.ident.[1]>highconfident)){output <- 1}
      else(output <- 0.5)
    }
    else(output <- 0)
  }
  else(output <- 0)
  return(output)
}

