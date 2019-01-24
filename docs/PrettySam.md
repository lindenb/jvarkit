# PrettySam

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Pretty SAM alignments


## Usage

```
Usage: prettysam [options] Files
  Options:
    -cN, --collapse-N
      collapse cigar operator 'N'
      Default: false
    -cD, --collapse-ND
      collapse cigar operator 'D'
      Default: false
    -color, --colors
      using ansi escape colors
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -kg, --knowngenes
      [20171219]UCSC knownGene File/URL. The knowGene format is a compact 
      alternative to GFF/GTF because one transcript is described using only 
      one line.	Beware chromosome names are formatted the same as your 
      REFERENCE. A typical KnownGene file is 
      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz 
      .If you only have a gff file, you can try to generate a knownGene file 
      with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html) 
      Memory intensive: the file is loaded in memory.
    -nA, --no-alignment
      hide alignment
      Default: false
    -nT, --no-attributes
      hide attributes table
      Default: false
    -nC, --no-cigar
      hide cigar string
      Default: false
    -noclip, --noclip, --no-clipping
      hide clipped bases
      Default: false
    -nP, --no-program-record
      hide program record
      Default: false
    -nR, --no-read-group
      hide read group
      Default: false
    -nS, --no-suppl
      hide supplementary alignements
      Default: false
    --no-unicode
      disable unicode to display ascii histogram
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -r, -R, --reference
      Indexed Genome Reference. It can be a the path to fasta file that must 
      be indexed with samtools faidx and with picard CreateSequenceDictionary. 
      It can also be a BioDAS dsn url like 
      `http://genome.cse.ucsc.edu/cgi-bin/das/hg19/` . BioDAS references are 
      slower, but allow to work without a local reference file.
    --trim
      trim long string to this length. <1 = do not trim.
      Default: 50
    -u, --unstranslated
      [20171219]Show untranslated regions (used with option -kg)
      Default: false
    -V, --variant, --vcf
      [20171220]Show VCF data. VCf must be indexed. if VCF has no genotype, 
      the variant positions are shown, otherwise, the genotypes associated to 
      the read will be show.
    --version
      print version and exit

```


## Keywords

 * sam
 * bam


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew prettysam
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/PrettySam.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/sam2tsv/PrettySam.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/PrettySamTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/sam2tsv/PrettySamTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **prettysam** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * PrettySam : a SAM/BAM prettifier. Lindenbaum & al. 2018. figshare. [https://doi.org/10.6084/m9.figshare.5853798.v1](https://doi.org/10.6084/m9.figshare.5853798.v1)


## About long reads

The htsjdk currently (2012-12-16) doesn't support more than 65635 cigar operations.

## Example

```

$ java -jar dist/prettysam.jar -R ref.fa S1.bam

>>>>> 1
          Read-Name : rotavirus_1_317_5:0:0_7:0:0_2de
               Flag : 99
             read paired : 1
             proper pair : 2
     mate reverse strand : 32
           first of pair : 64
               MAPQ : 60
             Contig : rotavirus  (index:0)
              Start : 1
                End : 70
             Strand : +
        Insert-Size : 317
        Mate-Contig : rotavirus  (index:0)
         Mate-Start : 248
        Mate-Strand : -
         Read Group : 
                      ID : S1
                      SM : S1
        Read-Length : 70
              Cigar : 70M (N=1)
           Sequence : 
                Read (0) : GGCTTTTAAT GCTTTTCAGT GGTTGCTGCT CAATATGGCG TCAACTCAGC AGATGGTCAG
                     Mid : |||||||||| |||||||||| |||||||||| ||| |||| | || ||||||| ||||||| ||
                 Ref (1) : GGCTTTTAAT GCTTTTCAGT GGTTGCTGCT CAAGATGGAG TCTACTCAGC AGATGGTAAG
                      Op : MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM
                    Qual : ++++++++++ ++++++++++ ++++++++++ ++++++++++ ++++++++++ ++++++++++
                 Ref-Pos : 1          11         21         31         41         51        

               Read (60) : CTCTAATATT
                     Mid : ||||| ||||
                Ref (61) : CTCTATTATT
                      Op : MMMMMMMMMM
                    Qual : ++++++++++
                 Ref-Pos : 61        

               Tags : 
                      MD :  33G4A3T14A7T4   "String for mismatching positions"
                      NM :  5   "Edit distance to the reference"
                      AS :  45   "Alignment score generated by aligner"
                      XS :  0   "Reserved for end users"
<<<<< 1

>>>>> 2
          Read-Name : rotavirus_1_535_4:0:0_4:0:0_1a6
               Flag : 163
             read paired : 1
             proper pair : 2
     mate reverse strand : 32
          second of pair : 128
               MAPQ : 60
             Contig : rotavirus  (index:0)
              Start : 1
                End : 70
             Strand : +
        Insert-Size : 535
        Mate-Contig : rotavirus  (index:0)
         Mate-Start : 466
        Mate-Strand : -
         Read Group : 
                      ID : S1
                      SM : S1
        Read-Length : 70
              Cigar : 70M (N=1)
           Sequence : 
                Read (0) : GGCTTTTACT GCTTTTCAGT GGTTGCTTCT CAAGATGGAG TGTACTCATC AGATGGTAAG
                     Mid : |||||||| | |||||||||| ||||||| || |||||||||| | |||||| | ||||||||||
                 Ref (1) : GGCTTTTAAT GCTTTTCAGT GGTTGCTGCT CAAGATGGAG TCTACTCAGC AGATGGTAAG
                      Op : MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM
                    Qual : ++++++++++ ++++++++++ ++++++++++ ++++++++++ ++++++++++ ++++++++++
                 Ref-Pos : 1          11         21         31         41         51        

               Read (60) : CTCTATTATT
                     Mid : ||||||||||
                Ref (61) : CTCTATTATT
                      Op : MMMMMMMMMM
                    Qual : ++++++++++
                 Ref-Pos : 61        

               Tags : 
                      MD :  8A18G13C6G21   "String for mismatching positions"
                      NM :  4   "Edit distance to the reference"
                      AS :  50   "Alignment score generated by aligner"
                      XS :  0   "Reserved for end users"
<<<<< 2
```

## Displaying genes:

```
$ wget -O knownGene.txt.gz "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz"
$ wget -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esG1bAlnRep1.bam" |\
	java -jar dist/prettysam.jar  -R "http://genome.cse.ucsc.edu/cgi-bin/das/hg19/" -kg knownGene.txt.gz

(...)

>>>>> 724
          Read-Name : SOLEXA-1GA-2_2_FC20EMB:5:4:549:957
               Flag : 0
               MAPQ : 25
             Contig : chr1  (index:0)
              Start : 1,115,507
                End : 1,115,542
             Strand : -->
        Read-Length : 36
              Cigar : 36M (N=1)
           Sequence : 
                Read (0) : GGCCGGACCT GGAGGGGGCA GAAAGAGCCT CCGCAA
                  Middle : |||||||||| |||||||||| |||||||||| | || |
         Ref (1,115,507) : ggccggacct ggagggggca gaaagagcct ctgcca
          Cigar-Operator : MMMMMMMMMM MMMMMMMMMM MMMMMMMMMM MMMMMM
                    Qual : TRGI^WDHH` RHIMOILGEF G@C>D>GC@= D@=A>@
                 Ref-Pos : 1115507    1115517    1115527    1115537   
      uc001acy.2(+) type : EEEEEEEEEE EEEEEEEEEE EEEEEEEEEE EEEEEE
                cDNA-Pos : 293        303        313        323       
                 Pep-Pos : 98         101        105        108       
              Amino-Acid : lyProAspLe uGluGlyAla GluArgAlaS erAlaT
      uc010nyg.1(+) type : EEEEEEEEEE EEEEEEEEEE EEEEEEEEEE EEEEEE
                cDNA-Pos : 293        303        313        323       
                 Pep-Pos : 98         101        105        108       
              Amino-Acid : lyProAspLe uGluGlyAla GluArgAlaS erAlaT
      uc001acz.2(+) type : EEEEEEEEEE EEEEEEEEEE EEEEEEEEEE EEEEEE
                cDNA-Pos : 74         84         94         104       
                 Pep-Pos : 25         28         32         35        
              Amino-Acid : lyProAspLe uGluGlyAla GluArgAlaS erAlaT

               Tags : 
                      X1 :           1   "Reserved for end users"
                      MD :      31C2A1   "String for mismatching positions"
                      NM :           1   "Edit distance to the reference"
<<<<< 724
```


## Screenshots

https://twitter.com/yokofakun/status/943132603539914752

![https://twitter.com/yokofakun/status/943132603539914752](https://pbs.twimg.com/media/DRatlqCWkAACABt.jpg)

https://twitter.com/yokofakun/status/942688906620887040

![https://twitter.com/yokofakun/status/942688906620887040](https://pbs.twimg.com/media/DRUaf8jWkAAArio.jpg)

https://twitter.com/yokofakun/status/941775073156968451

![https://twitter.com/yokofakun/status/941775073156968451](https://pbs.twimg.com/media/DRHbAN9XUAEFqdg.jpg)


