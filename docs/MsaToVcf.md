# MsaToVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Getting a VCF file from a CLUSTAW or a FASTA alignment. 


## DEPRECATED

use https://github.com/sanger-pathogens/snp_sites

## Usage

```
Usage: msa2vcf [options] Files
  Options:
    -R, --REF
      reference name used for the CHROM column. Optional
      Default: chrUn
    -a, --allsites
      print all sites
      Default: false
    -c, --consensus
      ruse this sequence as CONSENSUS
    -f, --fasta
      save computed fasta sequence in this file.
    -m, --haploid
      haploid output
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --version
      print version and exit

```


## Keywords

 * vcf
 * snp
 * msa
 * alignment



## See also in Biostars

 * [https://www.biostars.org/p/94573](https://www.biostars.org/p/94573)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew msa2vcf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/msa2vcf/MsaToVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/msa2vcf/MsaToVcf.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **msa2vcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



Deprecated: use https://github.com/sanger-pathogens/snp_sites , though some people told me they still use it for misc reasons.


## Motivation

Getting a VCF file from a CLUSTAW alignment. See also http://www.biostars.org/p/94573/

input is a clustalw file like: https://github.com/biopython/biopython/blob/master/Tests/Clustalw/opuntia.aln


## Cited-in


  * 'Differential distribution of Neandertal genomic signatures in human mitochondrial haplogroups'. 2017. Renata C Ferreira, Camila R Rodrigues, James R Broach, View ORCID ProfileMarcelo RS Briones. doi: [https://doi.org/10.1101/190363]([https://doi.org/10.1101/190363)
  * 'Pleiotropic effects of regulatory variation in tan result in correlation of two pigmentation traits in Drosophila melanogaster'. 2018. Molecular Ecology. Lukas Endler, Jean‐Michel Gibert, Viola Nolte, Christian Schlotterer. doi: 10.1111/mec.14781
  * 'Two key events associated with a transposable element burst occurred during rice domestication'  https://doi.org/10.1101/405290  
  * 'Tracking the origin of two genetic components associated with transposable element bursts in domesticated rice'. Nature Communicationsvolume 10, Article number: 641 (2019) https://www.nature.com/articles/s41467-019-08451-3

## Example

```bash
$ curl https://raw.github.com/biopython/biopython/master/Tests/Clustalw/opuntia.aln

CLUSTAL W (1.81) multiple sequence alignment


gi|6273285|gb|AF191659.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273284|gb|AF191658.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273287|gb|AF191661.1|AF191      TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273286|gb|AF191660.1|AF191      TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273290|gb|AF191664.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273289|gb|AF191663.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
gi|6273291|gb|AF191665.1|AF191      TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAA
                                    ******* **** *************************************

gi|6273285|gb|AF191659.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273284|gb|AF191658.1|AF191      TATATATA--------ATATATTTCAAATTTCCTTATATACCCAAATATA
gi|6273287|gb|AF191661.1|AF191      TATATA----------ATATATTTCAAATTTCCTTATATATCCAAATATA
gi|6273286|gb|AF191660.1|AF191      TATATA----------ATATATTTATAATTTCCTTATATATCCAAATATA
gi|6273290|gb|AF191664.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273289|gb|AF191663.1|AF191      TATATATATA------ATATATTTCAAATTCCCTTATATATCCAAATATA
gi|6273291|gb|AF191665.1|AF191      TATATATATATATATAATATATTTCAAATTCCCTTATATATCCAAATATA
                                    ******          ********  **** ********* *********

gi|6273285|gb|AF191659.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGT
gi|6273284|gb|AF191658.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273287|gb|AF191661.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273286|gb|AF191660.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273290|gb|AF191664.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
gi|6273289|gb|AF191663.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTAT
gi|6273291|gb|AF191665.1|AF191      AAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGT
                                    ************************************ *********** *

gi|6273285|gb|AF191659.1|AF191      ACCAGA
gi|6273284|gb|AF191658.1|AF191      ACCAGA
gi|6273287|gb|AF191661.1|AF191      ACCAGA
gi|6273286|gb|AF191660.1|AF191      ACCAGA
gi|6273290|gb|AF191664.1|AF191      ACCAGA
gi|6273289|gb|AF191663.1|AF191      ACCAGA
gi|6273291|gb|AF191665.1|AF191      ACCAGA
                                    ******
```
generate the VCF

```
$ curl https://raw.github.com/biopython/biopython/master/Tests/Clustalw/opuntia.aln" |\
  java -jar dist/msa2vcf.jar

##fileformat=VCFv4.1
##Biostar94573CmdLine=
##Biostar94573Version=ca765415946f3ed0827af0773128178bc6aa2f62
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth.">
##contig=<ID=chrUn,length=156>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	gi|6273284|gb|AF191658.1|AF191	gi|6273285|gb|AF191659.1|AF191	gi|6273286|gb|AF191660.1|AF191	gi|6273287|gb|AF191661.1|AF191	gi|6273289|gb|AF191663.1|AF191	gi|6273290|gb|AF191664.1|AF191	gi|6273291|gb|AF191665.1|AF191
chrUn	8	.	T	A	.	.	DP=7	GT:DP	0:1	0:1	1:1	0:1	0:1	0:1	0:1
chrUn	13	.	A	G	.	.	DP=7	GT:DP	0:1	0:1	0:1	0:1	1:1	1:1	1:1
chrUn	56	.	ATATATATATA	ATA,A,ATATA	.	.	DP=7	GT:DP	1:1	2:1	2:1	2:1	3:1	3:1	0:1
chrUn	74	.	TCA	TAT	.	.	DP=7	GT:DP	0:1	0:1	1:1	0:1	0:1	0:1	0:1
chrUn	81	.	T	C	.	.	DP=7	GT:DP	0:1	0:1	0:1	0:1	1:1	1:1	1:1
chrUn	91	.	T	C	.	.	DP=7	GT:DP	1:1	1:1	0:1	0:1	0:1	0:1	0:1
chrUn	137	.	T	C	.	.	DP=7	GT:DP	0:1	1:1	0:1	0:1	0:1	0:1	0:1
chrUn	149	.	G	A	.	.	DP=7	GT:DP	0:1	0:1	0:1	0:1	1:1	0:1	0:1
```

