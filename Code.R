library("seqinr")


leprae <- read.fasta(file = "Q9CD83.fasta")
ulcerans <- read.fasta(file = "Q9CD85.fasta")

lepraeseq <- leprae[[1]]
ulceransseq <- ulcerans[[1]]

lepraeseq # Display the contents of the vector "lepraeseq"

#Comparing two sequences using a dotplot
dotPlot(lepraeseq, ulceransseq)

#Pairwise global alignment of DNA sequences using the Needleman-Wunsch algorithm

library(Biostrings)

sigma <- nucleotideSubstitutionMatrix(match = 2, mismatch = -1, baseOnly = TRUE)
sigma # Print out the matrix

s1 <- "GAATTC"
s2 <- "GATTA"

globalAligns1s2 <- pairwiseAlignment(s1, s2, substitutionMatrix = sigma, gapOpening = -2,
gapExtension = -8, scoreOnly = FALSE)

globalAligns1s2

#Pairwise global alignment of protein sequences using the Needleman-Wunsch algorithm

data(BLOSUM50)
BLOSUM50 # Print out the data

#You can get a list of the available scoring matrices that come with the Biostrings package by using the data() function, which takes as an argument the name of the package for which you want to know the data sets that come with it:
data(package="Biostrings")

#To find the optimal global alignment between the protein sequences “PAWHEAE” and “HEAGAWGHEE” using the Needleman-Wunsch algorithm using the BLOSUM50 matrix, we type:


data(BLOSUM50)
s3 <- "PAWHEAE"
s4 <- "HEAGAWGHEE"
globalAligns3s4 <- pairwiseAlignment(s3, s4, substitutionMatrix = "BLOSUM50", gapOpening = -2,
gapExtension = -8, scoreOnly = FALSE)
globalAligns3s4 # Print out the optimal global alignment and its score

#Aligning UniProt sequences

lepraeseqstring <- c2s(lepraeseq)     # Make a string that contains the sequence in "lepraeseq"
ulceransseqstring <- c2s(ulceransseq) # Make a string that contains the sequence in "ulceransseq"


#Furthermore, pairwiseAlignment() requires that the sequences be stored as uppercase characters. Therefore, if they are not already in uppercase, we need to use the toupper() function to convert lepraeseqstring and ulceransseqstring to uppercase:
lepraeseqstring <- toupper(lepraeseqstring)
ulceransseqstring <- toupper(ulceransseqstring)
lepraeseqstring # Print out the content of "lepraeseqstring"

#We can now align the the M. leprae and M. ulcerans chorismate lyase protein sequences using the pairwiseAlignment() function:



globalAlignLepraeUlcerans <- pairwiseAlignment(lepraeseqstring, ulceransseqstring,
  substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
globalAlignLepraeUlcerans # Print out the optimal global alignment and its score


#Viewing a long pairwise alignment

printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
  {
     require(Biostrings)           # This function requires the Biostrings package
     seq1aln <- pattern(alignment) # Get the alignment for the first sequence
     seq2aln <- subject(alignment) # Get the alignment for the second sequence
     alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
     starts  <- seq(1, alnlen, by=chunksize)
     n       <- length(starts)
     seq1alnresidues <- 0
     seq2alnresidues <- 0
     for (i in 1:n) {
        chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
        chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
        # Find out how many gaps there are in chunkseq1aln:
        gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
        # Find out how many gaps there are in chunkseq2aln:
        gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
        # Calculate how many residues of the first sequence we have printed so far in the alignment:
        seq1alnresidues <- seq1alnresidues + chunksize - gaps1
        # Calculate how many residues of the second sequence we have printed so far in the alignment:
        seq2alnresidues <- seq2alnresidues + chunksize - gaps2
        if (returnlist == 'FALSE')
        {
           print(paste(chunkseq1aln,seq1alnresidues))
           print(paste(chunkseq2aln,seq2alnresidues))
           print(paste(' '))
        }
     }
     if (returnlist == 'TRUE')
     {
        vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
        vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
        mylist <- list(vector1, vector2)
        return(mylist)
     }
}


printPairwiseAlignment(globalAlignLepraeUlcerans, 60)

#Pairwise local alignment of protein sequences using the Smith-Waterman algorithm

localAlignLepraeUlcerans <- pairwiseAlignment(lepraeseqstring, ulceransseqstring,
  substitutionMatrix = BLOSUM50, gapOpening = -2, gapExtension = -8, scoreOnly = FALSE, type="local")

localAlignLepraeUlcerans # Print out the optimal local alignment and its score

printPairwiseAlignment(localAlignLepraeUlcerans, 60)

#Calculating the statistical significance of a pairwise global alignment
#We have seen that when we align the ‘PAWHEAE’ and ‘HEAGAWGHEE’ protein sequences, they


generateSeqsWithMultinomialModel <- function(inputsequence, X)
  {
     # Change the input sequence into a vector of letters
     require("seqinr") # This function requires the SeqinR package.
     inputsequencevector <- s2c(inputsequence)
     # Find the frequencies of the letters in the input sequence "inputsequencevector":
     mylength <- length(inputsequencevector)
     mytable <- table(inputsequencevector)
     # Find the names of the letters in the sequence
     letters <- rownames(mytable)
     numletters <- length(letters)
     probabilities <- numeric() # Make a vector to store the probabilities of letters
     for (i in 1:numletters)
     {
        letter <- letters[i]
        count <- mytable[[i]]
        probabilities[i] <- count/mylength
     }
     # Make X random sequences using the multinomial model with probabilities "probabilities"
     seqs <- numeric(X)
     for (j in 1:X)
     {
        seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
        seq <- c2s(seq)
        seqs[j] <- seq
     }
     # Return the vector of random sequences
     return(seqs)
  }


randomseqs <- generateSeqsWithMultinomialModel('PAWHEAE',1000)
randomseqs[1:10] # Print out the first 10 random sequences

s4 <- "HEAGAWGHEE"

pairwiseAlignment(s4, randomseqs[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
  gapExtension = -8, scoreOnly = FALSE)

pairwiseAlignment(s4, randomseqs[1], substitutionMatrix = "BLOSUM50", gapOpening = -2,
  gapExtension = -8, scoreOnly = TRUE)

randomscores <- double(1000) # Create a numeric vector with 1000 elements
for (i in 1:1000)
  {
     score <- pairwiseAlignment(s4, randomseqs[i], substitutionMatrix = "BLOSUM50",
       gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)
     randomscores[i] <- score
  }


hist(randomscores, col="red") # Draw a red histogram

sum(randomscores >= -5)






