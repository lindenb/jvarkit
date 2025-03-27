# Sam2Tsv

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Prints the SAM alignments as a TAB delimited file.


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar sam2tsv  [options] Files

Usage: sam2tsv [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --left-align
      Left-aligns any indels in the read data contained in a BAM or CRAM file. 
      The same indel can often be placed at multiple positions and still 
      represent the same haplotype. While it is a commonly used convention to 
      place an indel at the left-most position
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Indexed fasta Reference file. This file must be indexed with samtools 
      faidx and with picard/gatk CreateSequenceDictionary or samtools dict
    --regions
      Limit analysis to this interval. A source of intervals. The following 
      suffixes are recognized: vcf, vcf.gz bed, bed.gz, gtf, gff, gff.gz, 
      gtf.gz.Otherwise it could be an empty string (no interval) or a list of 
      plain interval separated by '[ \t\n;,]'
    -N, --skip-N
      Skip 'N' operator
      Default: false
    --validation-stringency
      SAM Reader Validation Stringency
      Default: LENIENT
      Possible Values: [STRICT, LENIENT, SILENT]
    --version
      print version and exit

```


## Keywords

 * sam
 * bam
 * table
 * cram
 * tsv



## See also in Biostars

 * [https://www.biostars.org/p/157232](https://www.biostars.org/p/157232)
 * [https://www.biostars.org/p/59647](https://www.biostars.org/p/59647)
 * [https://www.biostars.org/p/253828](https://www.biostars.org/p/253828)
 * [https://www.biostars.org/p/264875](https://www.biostars.org/p/264875)
 * [https://www.biostars.org/p/277493](https://www.biostars.org/p/277493)



## NF-CORE

![nfcorelogo](https://avatars.githubusercontent.com/u/35520196?s=32&v=4) This program is available in nf-core at [https://nf-co.re/modules/jvarkit_sam2tsv.html](https://nf-co.re/modules/jvarkit_sam2tsv.html)


## Creation Date

20170712

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2Tsv.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2TsvTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/Sam2TsvTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam2tsv** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example
 


```
$  java -jar dist/sam2tsv.jar -R src/test/resources/toy.fa src/test/resources/toy.bam 


#Read-Name  Flag  MAPQ  CHROM  READ-POS0  READ-BASE  READ-QUAL  REF-POS1  REF-BASE  CIGAR-OP
r001        163   30    ref    0          T          .          7         T         M
r001        163   30    ref    1          T          .          8         T         M
r001        163   30    ref    2          A          .          9         A         M
r001        163   30    ref    3          G          .          10        G         M
r001        163   30    ref    4          A          .          11        A         M
r001        163   30    ref    5          T          .          12        T         M
r001        163   30    ref    6          A          .          13        A         M
r001        163   30    ref    7          A          .          14        A         M
r001        163   30    ref    8          A          .          .         .         I
r001        163   30    ref    9          G          .          .         .         I
r001        163   30    ref    10         A          .          .         .         I
r001        163   30    ref    11         G          .          .         .         I
r001        163   30    ref    12         G          .          15        G         M
r001        163   30    ref    13         A          .          16        A         M
r001        163   30    ref    14         T          .          17        T         M
r001        163   30    ref    15         A          .          18        A         M
r001        163   30    ref    .          .          .          19        G         D
r001        163   30    ref    16         C          .          20        C         M
r001        163   30    ref    17         T          .          21        T         M
r001        163   30    ref    18         G          .          22        G         M
r002        0     30    ref    0          A          .          8         T         S
r002        0     30    ref    1          A          .          .         .         I
r002        0     30    ref    2          A          .          .         .         I
r002        0     30    ref    3          A          .          9         A         M
r002        0     30    ref    4          G          .          10        G         M
r002        0     30    ref    5          A          .          11        A         M
r002        0     30    ref    6          T          .          12        T         M
r002        0     30    ref    7          A          .          13        A         M
r002        0     30    ref    8          A          .          14        A         M
r002        0     30    ref    .          .          .          .         .         P
r002        0     30    ref    9          G          .          .         .         I
r002        0     30    ref    .          .          .          .         .         P
r002        0     30    ref    10         G          .          .         .         I
r002        0     30    ref    11         G          .          15        G         M
r002        0     30    ref    12         A          .          16        A         M
(...)

```


## Example 2

sam2tsv can read data from a linux pipe.

```
samtools view -h input.bam | java -jar dist/sam2tsv.jar
```



### Citations


Sam2tsv was cited in :

  * "Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements" . McCoy RC, Taylor RW, Blauwkamp TA, Kelley JL, Kertesz M, et al. (2014) Illumina TruSeq Synthetic Long-Reads Empower De Novo Assembly and Resolve Complex, Highly-Repetitive Transposable Elements. PLoS ONE 9(9): e106689. doi: 10.1371/journal.pone.0106689  http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0106689
  * "High-Throughput Identification of Genetic Variation Impact on pre-mRNA Splicing Efficiency". Scott I Adamson, Lijun Zhan, Brenton R Graveley. doi: [https://doi.org/10.1101/191122](https://doi.org/10.1101/191122).
  * "Linkage of A-to-I RNA editing in metazoans and the impact on genome evolution "  Molecular Biology and Evolution, msx274, https://doi.org/10.1093/molbev/msx274
  * "Vex-seq: high-throughput identification of the impact of genetic variation on pre-mRNA splicing efficiency" Genome Biology201819:71 https://doi.org/10.1186/s13059-018-1437-x
  * "Accurate detection of m6A RNA modifications in native RNA sequences" Huanle Liu, Oguzhan Begik, Morghan C Lucas, Christopher E Mason, Schraga Schwartz, John S Mattick, Martin A Smith, Eva Maria Novoa bioRxiv 525741; doi: https://doi.org/10.1101/525741 
  * "DART-seq: an antibody-free method for global m6A detection"  Nature Methods https://doi.org/10.1038/s41592-019-0570-0
  * "Thiouridine-to-Cytidine Conversion Sequencing (TUC-Seq) to Measure mRNA Transcription and Degradation Rates" The Eukaryotic RNA Exosome. Nov 2019. https://doi.org/10.1007/978-1-4939-9822-7_10
  * "Evolutionary forces on A-to-I RNA editing revealed by sequencing individual honeybee drones". Yuange Duan, Shengqian Dou, Jiaxing Huang, Eli Eisenberg, Jian Lu . 2020 . https://doi.org/10.1101/2020.01.15.907287
  * "Sci-fate characterizes the dynamics of gene expression in single cells" (2020) Nat Biotechnol (2020). https://doi.org/10.1038/s41587-020-0480-9
  * "Mutations in virus-derived small RNAs" Nigam Deepti; LaTourrette, Katherine; Garcia-Ruiz, Hernan. Scientific Reports (Nature Publisher Group); London Vol. 10, Iss. 1,  (2020). DOI:10.1038/s41598-020-66374-2 
  * Qiu, Q., Hu, P., Qiu, X. et al. Massively parallel and time-resolved RNA sequencing in single cells with scNT-seq. Nat Methods (2020). https://doi.org/10.1038/s41592-020-0935-4
  * FUJIKURA, K. et al. Multiregion whole-exome sequencing of intraductal papillary mucinous neoplasms reveals frequent somatic KLF4 mutations predominantly in low-grade regions. Gut, [s. l.], 2020. DOI 10.1136/gutjnl-2020-321217
  * Gao, Y., Liu, X., Wu, B. et al. Quantitative profiling of N6-methyladenosine at single-base resolution in stem-differentiating xylem of Populus trichocarpa using Nanopore direct RNA sequencing. Genome Biol 22, 22 (2021). https://doi.org/10.1186/s13059-020-02241-7
  * Liu H., Begik O., Novoa E.M. (2021) EpiNano: Detection of m6A RNA Modifications Using Oxford Nanopore Direct RNA Sequencing. In: McMahon M. (eds) RNA Modifications. Methods in Molecular Biology, vol 2298. Humana, New York, NY. https://doi.org/10.1007/978-1-0716-1374-0_3
  * Yang & al. "Sequencing 5-Formyluracil in Genomic DNA at Single-Base Resolution"  (2021) Analytical Chemistry doi: 10.1021/acs.analchem.1c03339
  * Xiao Shu, Chenyang Huang, Tengwei Li, Jie Cao, Jianzhao Liu,.a6A-seq: N6-allyladenosine-based cellular messenger RNA metabolic labelling and sequencing,.Fundamental Research,.2023,.,.ISSN 2667-3258, https://doi.org/10.1016/j.fmre.2023.04.010..(https://www.sciencedirect.com/science/article/pii/S2667325823001267).
  * Xu, Z., Sziraki, A., Lee, J. et al. Dissecting key regulators of transcriptome kinetics through scalable single-cell RNA profiling of pooled CRISPR screens. Nat Biotechnol (2023). https://doi.org/10.1038/s41587-023-01948-9
  * Braton & al (2023) HIVepsilon-seq-scalable characterization of intact persistent proviral HIV reservoirs in women. Virology. https://doi.org/10.1128/jvi.00705-23
  * Yin, K., Xu, Y., Guo, Y. et al. Dyna-vivo-seq unveils cellular RNA dynamics during acute kidney injury via in vivo metabolic RNA labeling-based scRNA-seq. Nat Commun 15, 9866 (2024). https://doi.org/10.1038/s41467-024-54202-4
  * Heng, J. J. W. (2024). Development of sequencing methods to illuminate the epitranscriptome. Doctoral thesis, Nanyang Technological University, Singapore. https://hdl.handle.net/10356/182139
  * Liu, C., Liang, H., Wan, A.H. et al. Decoding the m6A epitranscriptomic landscape for biotechnological applications using a direct RNA sequencing approach. Nat Commun 16, 798 (2025). https://doi.org/10.1038/s41467-025-56173-6


