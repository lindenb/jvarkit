---
title: 'BamWithoutBai: title'
tags:
  - java
  - hts
  - bioinformatics
  - bai
  - bam
  - bgzf
authors:
  - name: Pierre Lindenbaum
    orcid: 0000-0003-0148-9787
    affiliation: "1" 
  - name: Tom White
    affiliation: 2
  - name: Richard Redon
    orcid: 0000-0001-7751-2280
    affiliation: "1" 
affiliations:
 - name: L'institut du thorax, INSERM, CNRS, Université de Nantes, F-44000 Nantes, France.
   index: 1
 - name: undefined
   index: 2
date: 19 December 2019
bibliography: ../paper.bib


---

# Summary

*BAM* (*Binary Alignment/Map format* ) is a [standard](https://samtools.github.io/hts-specs/SAMv1.pdf) binary format for storing hight-troughtpout sequencing data [@pmid19505943].


*BAM* files are compressed using the *BGZF* format. The **BGZF** format is a compression format implemented on top of the standard gzip file format while allowing efficient random access to the BAM file for indexed queries. A Bam Index (*Bai*) achieves a fast retrieval of alignments overlapping a specified region without going through the whole alignment. It also works through the web. 


Encode

The quick brown fox jumps over the lazy dog [@htsjdk] , the quick brown fox jumps over the lazy dog [@pmid19505943] ,  the quick brown fox jumps over the lazy dog [@jvarkit].

```
$ cat interval.bed
chr22	41697506	41756151

$ time wget -q -O - "https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam" |\
	samtools view -L interval.bed -O BAM -owget.bam - 

real	17m48,710s
user	8m7,346s
sys	0m51,719s


$ samtools view wget.bam | sha1sum 
5fc261d051592edd791309211727ebbd0e3de909  -

time java -jar dist/bamwithoutbai.jar  -o nobai.bam -r "chr22:41697507-41756151" \
	 'https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam' 

real	0m24,848s
user	0m3,860s
sys	0m0,320s

$ samtools view nobai.bam | sha1sum 
5fc261d051592edd791309211727ebbd0e3de909  -
```



# Acknowledgements

We acknowledge contributions from Louis Bergelson (Broad Institute) and John Marshall during the genesis of this project.

# References

