setwd('C:/Users/lenda/Desktop/VUT_Ing/1_semestr/MPA-PRG/exercises/MPA_PRG_exercise_07')

library('Biostrings')

#######################################################################
#################### The Brute-Force Motif Search #####################
#######################################################################
# TASK 1
Score <- function(s,DNA,l){
  block <- DNAStringSet()
  for (seq_idx in 1:length(DNA)){
    block <- c(block, subseq(DNA[seq_idx],
                             start = s[seq_idx], width = l)) # extract subbseq
  }
  frequency <- consensusMatrix(block) # generate frequency matrix
  score <- 0

  for (seq_idx in 1:ncol(frequency)){
    score <- score + max(frequency[,seq_idx])
  }

  return(score)
}

# TASK 2
NextLeaf <- function(a, L, k){
  for (i in L:1){
    if (a[i] < k){
      a[i] <- a[i] + 1
      return(a)
    }
    a[i] <- 1
  }
  return(a)
}

# TASK 3
BFMotifSearch <- function(DNA, L, n, l){
  s <- rep(1, L)
  bestScore <- Score(s,DNA,l)
  while (TRUE){
    s <- NextLeaf(s, L, n-l+1)
    if (Score(s,DNA,l) > bestScore){
      bestScore <- Score(s,DNA,l)
      bestMotif <- s
    }
    if (all(s == rep(1, L))){
      return(bestMotif)
    }
  }
}

# DNA <- readDNAStringSet('seq_motif.fasta')
# s <- c(1,1,1,1,1)
# l <- 8
# Score(s, DNA, l)
# L <- length(DNA)
# n <- length(DNA[[1]])
# k <- n - l + 1
# NextLeaf(s, L, k)
# BFMotifSearch(DNA, L, n, l)


#######################################################################
################## The Branch-and-Bound Motif Search ##################
#######################################################################
# TASK 4
NextVertex <- function(a, i, L, k){
  if (i < L){
    a[i+1] <- 1
    return(list(a, i+1))
  }else{
    for (j in L:1){
      if (a[j] < k){
        a[j] <- a[j] + 1
        return(list(a,j))
      }
    }
  }
  return(list(a,0))
}

# TASK 5
ByPass <- function(a, i, L, k){
  for (j in i:1){
    if (a[j] < k){
      a[j] <- a[j] + 1
      return(list(a,j))
    }
  }
  return(list(a,0))
}

# TASK 6
BBMotidSearch <- function(DNA, t, n, l){
  s <- rep(1, t)
  bestScore <- 0
  i <- 1
  while (i > 0){
    if (i < t){
      optimisticScore <- Score2(s, i, DNA, l) + (t - i) * l
      if (optimisticScore < bestScore){
       pom <- ByPass(s, i, t, n-l+1)
        s <- pom[[1]]
        i <- pom[[2]]
      }else{
        pom <- NextVertex(s, i, t, n-l+1)
        s <- pom[[1]]
        i <- pom[[2]]
      }
    }else{
      if (Score2(s, t, DNA, l) > bestScore){
        bestScore <- Score2(s, t, DNA, l)
        bestMotif <- s
      }
      pom <- NextVertex(s, i, t, n-l+1)
      s <- pom[[1]]
      i <- pom[[2]]
    }
  }
  return(bestMotif)
}

Score2 <- function(s,i,DNA,l){
  block <- DNAStringSet()
  for (seq_idx in 1:i){
    block <- c(block, subseq(DNA[seq_idx],
                             start = s[seq_idx], width = l)) # extract subbseq
  }
  frequency <- consensusMatrix(block) # generate frequency matrix
  score <- 0

  for (seq_idx in 1:ncol(frequency)){
    score <- score + max(frequency[,seq_idx])
  }

  return(score)
}

#######################################################################

# DNA <- readDNAStringSet('seq_motif.fasta')
DNA <- readDNAStringSet('seq_score.fasta')

# s <- c(1,1,1,1,1)           # L-mer starting indexes
t <- length(DNA)            # number od sequences
n <- length(DNA[[1]])       # length of sequences
l <- 6                      # length of motif
# k <- n-l+1

BFMotifSearch(DNA, t, n, l) # Brute-Force Motif Search
BBMotidSearch(DNA, t, n, l) # Branch-and-Bound Motif Search


