# 03/14/2016
# Ryan R. Cupo
# Strucural Analysis of Skd3

# clears objects from working environment & initializes new
rm(list=ls())
objects()

# sets wokring directory
setwd("C:\\Users\\Ryan Cupo\\Documents\\Graduate School\\2016_Spring\\GCB 535\\Shorter Lab Bioinformatics Project")

# loads file from path into R studio
# load in ExAc file
DataSet <- read.csv(file.choose(), header = TRUE)
# load in amino acid reference table
aminoTable <- read.csv(file.choose(), header = TRUE)

# subsetting with refining ("|" symbol is "or")
DataSet <- subset(DataSet, (DataSet$Annotation == "missense"))

# subset relevant columns
DataSet <- subset(DataSet, select = c(Consequence, Annotation, Flags, Allele.Count, Number.of.Homozygotes, Allele.Frequency))

# find length of DataSet
lenDataSet <- nrow(DataSet)

# find length of aminoTable
lenAminoTable <- nrow(aminoTable)

# set counting vars
i <- 1
j <- 1

# loop through every cell in first column and remove non-numeric values to create column for residue location
while (i <= lenDataSet)
{
  # create and sort column for residue location
  DataSet$Location[i] <- gsub("[^a-zA-Z0-9]", "", DataSet$Consequence[i])
  DataSet$Location[i] <- gsub("[A-z]", "", DataSet$Location[i])
  # create and sort column for residue identity
  strCat <- gsub("[^a-zA-Z0-9]", "", DataSet$Consequence[i])
  strCat <- substring(strCat, seq(1,nchar(strCat),1), seq(1,nchar(strCat),1))
  strCat <- paste(strCat[2],strCat[3],strCat[4],sep='')
  DataSet$Residue[i] <- strCat
  # loop to compare 3 letter amino acid residues with database and write 1 letter equivalent
  while (j <= lenAminoTable)
  {
    if (DataSet$Residue[i] == aminoTable$Abbreviation3[j])
    {
      print(aminoTable$Abbreviation1[j])
      print(DataSet$Residue[i])
      DataSet$Residue1[i] <- toString(aminoTable$Abbreviation1[j])
    }
    j <- j + 1
  }
  j <- 1
  i <- i + 1
}

# convert Location to numeric type and then sort by residue
attach(DataSet)
Location <- as.numeric(Location)
DataSet <- DataSet[order(Location, Consequence),]
detach(DataSet)

# reset counting var
i <- 1

# sum repeat residues
while (i <= lenDataSet)
{
  # tests if there are multiple missense mutations at same residue
  loci <- DataSet$Location[i]
  lociplus1 <- DataSet$Location[i+1]
  test <- identical(loci, lociplus1)
  # uses test values to either sum counts for duplicates (test == TRUE) or progress through data
  if (test == TRUE)
  {
    # sum counter columns
    DataSet$Allele.Count[i+1] <- DataSet$Allele.Count[i] + DataSet$Allele.Count[i+1]
    DataSet$Number.of.Homozygotes[i+1] <- DataSet$Number.of.Homozygotes[i] + DataSet$Number.of.Homozygotes[i+1]
    DataSet$Allele.Frequency[i+1] <- DataSet$Allele.Frequency[i] + DataSet$Allele.Frequency[i+1]
    # set test value to TRUE that there are duplicates (for later sorting)
    DataSet$dup[i] <- TRUE
    # add to count var
    i <- i + 1
    # other thoughts on how to remove duplicates
    ## remove NaNs?
    ## maybe use 'try' to do test
    ## function like 'splice' idk what is called
    ## I could also make a new table that stores the processed result
  }
  else if (test == FALSE)
  {
    # set test value to FALSE (unique sum hit for sorting)
    DataSet$dup[i] <- FALSE
    # add to count var
    i <- i + 1
  }
}

# sort only non duplicate rows (DataSet$dup == FALSE)
DataSet <- DataSet[!(DataSet$dup==TRUE),]
## d<-d[!(d$A=="B" & d$E==0),]

# remove 'dup' column
DataSet <- subset(DataSet, select = c(Consequence, Annotation, Flags, Allele.Count, Number.of.Homozygotes, Allele.Frequency, Location, Residue, Residue1))

# write table as csv
write.csv(DataSet, file = "RyanCupo_RScriptOutput.csv", row.names = TRUE)

# subset dataset to remove all residues out of range of 3d model
DataSetNBD2 <- DataSet[which(as.numeric(DataSet$Location)>=308 & as.numeric(DataSet$Location)<=663),]

# find nrow datasetnbd2
lenDataSetNBD2 <- nrow(DataSetNBD2)

# reset counting var to 1
i <- 1

# set DataSetNBD2$Location to numeric and sort by order
DataSetNBD2$Location <- as.numeric(DataSetNBD2$Location)

# set residue 0 as start of NBD2 from SWISS-MODEL
while (i <= lenDataSetNBD2)
{
  DataSetNBD2$Location[i] <- DataSetNBD2$Location[i] - 307
  i <- i + 1
}

# write NBD2 table as csv
write.csv(DataSetNBD2, file = "RyanCupo_RScriptOutputNBD2.csv", row.names = TRUE)

# filter by frequency (allele count > 10)
DataSetFreq <- DataSet[which(as.numeric(DataSet$Allele.Count)>10),]

# filter by homozygotes
DataSetHomo <- DataSet[which(as.numeric(DataSet$Number.of.Homozygotes)>0),]

# to interface with PyMol I need to have a text list of all the residues in each group connected as (1+2+3+7+11+14)
posList <- DataSet$Location

# writes SNP list to .txt in dir
sink('SNPForPyMol.txt')
cat(posList, sep = '+')
sink()

# load in disease causing mutants list (3-methylglutaconic aciduria)
MUTForPyMol <- read.csv(file.choose(), header = TRUE)

# writes MUT list to .txt in dir
sink('MUTForPyMol.txt')
cat(MUTForPyMol$X408, sep = '+')
sink()

# writes FRQ list to .txt in dir
sink('FRQForPyMol.txt')
cat(DataSetFreq$Location, sep = '+')
sink()

# writes HOM list to .txt in dir
sink('HOMForPyMol.txt')
cat(DataSetHomo$Location, sep = '+')
sink()

# now we are ready to interface with PyMol (hopefully)