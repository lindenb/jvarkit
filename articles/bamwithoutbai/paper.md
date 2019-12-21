---
title: 'BamWithoutBai: fast retrieval of alignments from remote BAM without BAI index.'
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
 - name: Unaffiliated
   index: 2
date: 19 December 2019
bibliography: ../paper.bib


---

# Summary

Remote files storing sequencing data (*BAM*) can be queried by genomic intervals if they are associated with an index file (*BAI*). If the *BAI* is missing, the whole *BAM* must be downloaded and filtered on the client side. We have developed a tool named '*bamwithoutbai*' which is able to quickly download a genomic region without this associated index.

# Description

*BAM* (*Binary Alignment/Map format* ) is a [standard](https://samtools.github.io/hts-specs/SAMv1.pdf) binary format for storing high-throughput sequencing data [@pmid19505943]. *BAM* files are compressed using the *BGZF* format, a compression format implemented on top of the standard gzip file format while allowing efficient random access to the BAM file for indexed queries. An associated BAM Index (*BAI*) achieves the fast retrieval of alignments overlapping a specified region without going through the whole alignment. This technology also works through HTTP, making it possible to fetch regions from remote BAM files.

Some sequencing projects, like ENCODE [@pmid29126249], provide *BAM* files but don't provide the associated *BAI* files, making it difficult to retrieve some specific regions of the BAM files: users need to download the whole *BAM* file before filtering the overlapping reads with, for example, `samtools`.

We wrote a program named '*bamwithoutbai*' that performs a binary search on a remote *BAM*: it sends a HTTP `Range-Request` at some point of the BAM, tries to find the next *BGZF* block and uncompress the stream until it finds a valid *BAM* record. This operation is repeated until the region of interest is exhausted.

In our test, downloading and filtering the region `chr3:38548061-38649667` for the ENCODE *BAM* '[ENCFF741DEO.bam](https://www.encodeproject.org/files/ENCFF741DEO/@@download/ENCFF741DEO.bam)' (12 gigabytes) with `samtools`  [@pmid19505943] takes about 17 minutes, whereas it only takes 24 seconds using our software.

`BamWithoutBai` is a Java program, it uses the 'java library for high throughput sequencing' [@htsjdk] as well as some code from Hadoop-BAM project [@pmid22302568]. It's part of the 'jvarkit' package [@jvarkit] and is available under the MIT License at: [http://lindenb.github.io/jvarkit/BamWithoutBai.html](http://lindenb.github.io/jvarkit).

# Acknowledgements

We acknowledge the useful suggestions from [Louis Bergelson](https://github.com/samtools/htsjdk/issues/1445#issuecomment-565599459) and [John Marshall](https://twitter.com/jomarnz/status/1205532441353560066) during the genesis of this project.

# References
