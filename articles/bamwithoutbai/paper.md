---
title: 'BamWithoutBai: fast retrieval of alignments from remote BAM without index.'
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

*BAM* (*Binary Alignment/Map format* ) is a [standard](https://samtools.github.io/hts-specs/SAMv1.pdf) binary format for storing hight-troughtpout sequencing data [@pmid19505943]. *BAM* files are compressed using the *BGZF* format. The **BGZF** format is a compression format implemented on top of the standard gzip file format while allowing efficient random access to the BAM file for indexed queries. A Bam Index (*Bai*) achieves a fast retrieval of alignments overlapping a specified region without going through the whole alignment. It also works through the web. 


Some sequencing project like Encode[@pmid29126249] provides bam files but don't provide the associated BAI files, making it difficult to retrieve some specific regions of the BAM files: a user needs to download the whole BAM file before filtering the overlaping reads with samtools. 

We wrote a progam named 'bamwithoutbai' that run a kind of binary search on a remote BAM: it performs a HTTP Range-Request at some point of the bam, tries to find the next *BGZF* block and uncompress the stream until it finds a valid BAM record. This operation is executed until the region of interest is narrowed.

In our test, downloading and filtering a BAM with samtools takes about 17 minutes while using our software takes 24 seconds.

BamWithoutBai is a java program, it uses the 'java library for high throughput sequencing'[@htsjdk] as well as some code from Hadoop-Bam project[@pmid22302568]. , it's part of the 'jvarkit' package[@jvarkit] and is available under the MIT License at: [http://lindenb.github.io/jvarkit/BamWithoutBai.html](http://lindenb.github.io/jvarkit).


# Acknowledgements

We acknowledge contributions from Louis Bergelson (Broad Institute) and John Marshall during the genesis of this project.

# References

