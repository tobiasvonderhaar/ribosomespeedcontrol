require(seqinr)

yeastORFs <- read.fasta('orf_coding.fasta', as.string = TRUE)
codons <- read.csv('codons.csv')

firstTen <- as.data.frame(mat.or.vec(length(yeastORFs), 4))
colnames(firstTen) <- c('ORF', 'firstTenCodons', 'DecodingTime', 'P')

#extract the first 30 nt for each gene
for(n in 1:length(yeastORFs)) {
  firstTen$ORF[n] <- getName(yeastORFs[[n]])
  firstTen$firstTenCodons[[n]] <- toupper(c2s(getSequence(yeastORFs[[n]])[1:30]))
}

#collate passage times for the first ten codons
for(m in 1:dim(firstTen)[1]) {
  n = 1
  while(n <30) {
    firstTen[m,"DecodingTime"] = firstTen[m,"DecodingTime"] + codons[codons$codon == substr(firstTen$firstTenCodons[m], start = n, stop = n + 2),"decoding.time"]
    n = n + 3
  }
}

print('Calculated actual decoding times, now processing decoding time percentiles')

##generate random sequences, evaluate their decoding times, compare actual decoding time

#load uniweight variable for the permutation analysis
load('uniweight.RData')

for(gene in 1:5917) {
  sampleno = 1
  this.samplevector <- mat.or.vec(1,1)
  sampling_sufficient = FALSE
  samplepos = 1
  sampleno = 1
  
  #generate a vector with decoding times of sequences with random codon composition. Sample at least 500 times, or until the fractional position is < 1/5*sampleno

  #sample at least 500 times
  while(sampleno < 500) {
      this.samplevector[sampleno] <- 0
      ntno = 1
      this.permutation <- paste(permutation(s2c(as.character(firstTen[gene,"firstTenCodons"])), modele = 'syncodon', ucoweight = uniweight), collapse = "")
      this.permutation <- toupper(this.permutation)
      while(ntno <30) {
        this.samplevector[sampleno] = this.samplevector[sampleno] + codons[codons$codon == substr(this.permutation, start = ntno, stop = ntno + 2),"decoding.time"]
        ntno = ntno + 3
      }
      sampleno = sampleno + 1
  }
  ordered.samplevector <- sort(this.samplevector)
  percentile = min(which(firstTen[gene,"DecodingTime"] < ordered.samplevector)) / sampleno
  
 #test whether the number of samples is sufficient to determine the rank reliably
  if(percentile > (1 / sampleno * 5) && percentile < (1 - (1/sampleno * 5))) {
    sampling_sufficient = TRUE
  }
  
  #if sampling is insufficient, add more samples in batches of 100 and re-test whether sampling is sufficient
  while(sampling_sufficient == FALSE) {  
    for(counter in 1:100) {
        this.samplevector[sampleno] <- 0
        ntno = 1
        this.permutation <- paste(permutation(s2c(as.character(firstTen[gene,"firstTenCodons"])), modele = 'syncodon', ucoweight = uniweight), collapse = "")
        this.permutation <- toupper(this.permutation)
        while(ntno <30) {
          this.samplevector[sampleno] = this.samplevector[sampleno] + codons[codons$codon == substr(this.permutation, start = ntno, stop = ntno + 2),"decoding.time"]
          ntno = ntno + 3
        }
        sampleno = sampleno + 1
      }
      ordered.samplevector <- sort(this.samplevector)
      percentile = min(which(firstTen[gene,"DecodingTime"] < ordered.samplevector)) / sampleno
      if(percentile > (1 / sampleno * 5) && percentile < (1 - (1/sampleno * 5))) {
        sampling_sufficient = TRUE
      }
      if(sampleno > 500000) {
        sampling_sufficient = TRUE
        percentile = 1/500000
      }
  }
  firstTen[gene,"DecodingTimePercentile"] <- percentile
  if(ceiling(gene/10) == gene/10) {print(c(gene, as.character(firstTen[gene,"ORF"])))}
}

write.csv(firstTen, 'firstTenNew.csv')
  



