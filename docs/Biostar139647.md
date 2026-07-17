# Biostar139647

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert alignment in Fasta/Clustal format to SAM/BAM file


## Usage


This program is now part of the main `jvarkit` tool. See [jvarkit](JvarkitCentral.md) for compiling.


```
Usage: java -jar dist/jvarkit.jar biostar139647  [options] Files

Usage: biostar139647 [options] Files
  Options:
    --bamcompression
      Compression Level. 0: no compression. 9: max compression;
      Default: 5
    -F, --fullalignment
      by default heading/traling gaps are ignored ( as in read to reference 
      alignment). Use this option to treat them as meaningfull indels as in 
      full genome alignment. Optional
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -X, --matchandmismatch
      when a reference is selected, identical residue match (=) and mismatches 
      (X) are output. Default is using general match (M). Optional
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -R, --refname
      reference name injected in the SAM/BAM IF not present in the alignment. 
      Note that not defining a reference will produce no Insertion.  Optional
      Default: chrUn
    -S, --refselection
      label of the reference, if present in the input alignment. Both 
      Insertion and Deletions variations can be computed. Optional
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```


## Keywords

 * msa
 * sam
 * bam
 * clustal



## See also in Biostars

 * [https://www.biostars.org/p/139647](https://www.biostars.org/p/139647)


## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar139647.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar139647.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar139647Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar139647Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar139647** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

 
## Example

```bash
$ curl -sL "https://raw.githubusercontent.com/suryasaha/Pred_cutoff/60a6f980c9940dfb6e381c5394918f27cb14564f/data/Xylella-RpoH.aln" |\
  java -jardist/biostar139647.jar

@HD	VN:1.4	SO:unsorted
@SQ	SN:chrUn	LN:42
@PG	ID:0	VN:3a0c4ccb05e7492382e00328ac60951f215d9400	CL:(empty)	PN:Biostar139647
1	0	chrUn	1	60	42M	*	0	0	CATACTTGGTCATCGGTCGTGTCCTTGAAAGTGACTTGTTAA	*
2	0	chrUn	1	60	42M	*	0	0	TCTCTGAACCCCCTTGAAACCCCTACACTCAGCCATATATGC	*
3	0	chrUn	1	60	42M	*	0	0	TACCTTCGGGTCCTTGAAAATAGCGTCGCCGTGCTTATCTGT	*
4	0	chrUn	1	60	5M2D35M	*	0	0	TTGACAGCCGCTTGAGCAGGCGTCGGTCATCCCCACATTC	*
5	0	chrUn	1	60	18M1D9M1D13M	*	0	0	ATGCCTGGGTGGCTTGAAAGCTGGCGGCTTGCCCACATAC	*
6	0	chrUn	1	60	20M1D21M	*	0	0	TCAGTTTTATCGCTTGATATTCACTGAGACTGGCCACACAT	*

```

```
$ curl -sL "https://raw.github.com/biopython/biopython/master/Tests/Clustalw/opuntia.aln" 2> /dev/null | java -jar dist-1.128/biostar139647.jar 
@HD	VN:1.4	SO:unsorted
@SQ	SN:chrUn	LN:156
@PG	ID:0	VN:3a0c4ccb05e7492382e00328ac60951f215d9400	CL:(empty)	PN:Biostar139647
gi|6273285|gb|AF191659.1|AF191	0	chrUn	1	60	56M10D90M	*	0	0	TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCA
AATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCCATTGATTTAGTGTACCAGA*
gi|6273284|gb|AF191658.1|AF191	0	chrUn	1	60	58M8D90M	*	0	0	TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATAATATATTT
CAAATTTCCTTATATACCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA	*
gi|6273289|gb|AF191663.1|AF191	0	chrUn	1	60	60M6D90M	*	0	0	TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATAT
TTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTATACCAGA	*
gi|6273291|gb|AF191665.1|AF191	0	chrUn	1	60	156M	*	0	0	TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATATATATAATATATTT
CAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA	*
gi|6273287|gb|AF191661.1|AF191	0	chrUn	1	60	56M10D90M	*	0	0	TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTCA
AATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA*
gi|6273286|gb|AF191660.1|AF191	0	chrUn	1	60	56M10D90M	*	0	0	TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATAATATATTTAT
AATTTCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA*
gi|6273290|gb|AF191664.1|AF191	0	chrUn	1	60	60M6D90M	*	0	0	TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAGAAAGAATATATATATAATATAT
TTCAAATTCCCTTATATATCCAAATATAAAAATATCTAATAAATTAGATGAATATCAAAGAATCTATTGATTTAGTGTACCAGA	*
```


