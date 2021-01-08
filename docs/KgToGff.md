# KgToGff

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Convert UCSC knowGene file to gff3


## Usage

```
Usage: kg2gff [options] Files
  Options:
    --coding
      select coding transcript only.
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    --source
      gff source field
      Default: ucsc
    --version
      print version and exit

```


## Keywords

 * gff
 * gff3
 * ucsc


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew kg2gff
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210106

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/KgToGff.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/gtf/KgToGff.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gtf/KgToGffTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/gtf/KgToGffTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **kg2gff** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


# Warning

ouput is escaped for UTF8, some characters like ':' might be converted to a hexadecimal encoding.

# Example

```
$ head -n1  ~/jeter.kg | java -jar dist/kg2gff.jar   | sed 's/%3A/:/g'  

##gff-version 3.1.25
chr22	ucsc	gene	18317101	18507325	.	-	.	ID=gene:uc002znh.2.g1;Name=uc002znh.2.1;biotype=protein_coding;gene_id=uc002znh.2.g1
chr22	ucsc	mRNA	18317101	18507325	.	-	.	ID=transcript:uc002znh.2.g1;Parent=gene:uc002znh.2.g1;Name=uc002znh.2.1;biotype=protein_coding;transcript_id=uc002znh.2.t1
chr22	ucsc	exon	18317101	18317315	.	-	.	ID=uc002znh.2:E0;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E0
chr22	ucsc	CDS	18317216	18317315	.	-	1	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS0;protein_id=uc002znh.2
chr22	ucsc	exon	18324588	18324783	.	-	.	ID=uc002znh.2:E1;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E1
chr22	ucsc	CDS	18324588	18324783	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS1;protein_id=uc002znh.2
chr22	ucsc	exon	18347665	18347752	.	-	.	ID=uc002znh.2:E2;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E2
chr22	ucsc	CDS	18347665	18347752	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS2;protein_id=uc002znh.2
chr22	ucsc	exon	18348690	18348778	.	-	.	ID=uc002znh.2:E3;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E3
chr22	ucsc	CDS	18348690	18348778	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS3;protein_id=uc002znh.2
chr22	ucsc	exon	18354603	18354789	.	-	.	ID=uc002znh.2:E4;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E4
chr22	ucsc	CDS	18354603	18354789	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS4;protein_id=uc002znh.2
chr22	ucsc	exon	18368644	18368817	.	-	.	ID=uc002znh.2:E5;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E5
chr22	ucsc	CDS	18368644	18368817	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS5;protein_id=uc002znh.2
chr22	ucsc	exon	18369936	18369998	.	-	.	ID=uc002znh.2:E6;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E6
chr22	ucsc	CDS	18369936	18369998	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS6;protein_id=uc002znh.2
chr22	ucsc	exon	18370089	18370201	.	-	.	ID=uc002znh.2:E7;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E7
chr22	ucsc	CDS	18370089	18370201	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS7;protein_id=uc002znh.2
chr22	ucsc	exon	18371800	18371996	.	-	.	ID=uc002znh.2:E8;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E8
chr22	ucsc	CDS	18371800	18371996	.	-	1	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS8;protein_id=uc002znh.2
chr22	ucsc	exon	18374251	18374398	.	-	.	ID=uc002znh.2:E9;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E9
chr22	ucsc	CDS	18374251	18374398	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS9;protein_id=uc002znh.2
chr22	ucsc	exon	18376574	18376670	.	-	.	ID=uc002znh.2:E10;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E10
chr22	ucsc	CDS	18376574	18376670	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS10;protein_id=uc002znh.2
chr22	ucsc	exon	18378050	18378176	.	-	.	ID=uc002znh.2:E11;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E11
chr22	ucsc	CDS	18378050	18378176	.	-	1	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS11;protein_id=uc002znh.2
chr22	ucsc	exon	18379012	18379127	.	-	.	ID=uc002znh.2:E12;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E12
chr22	ucsc	CDS	18379012	18379127	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS12;protein_id=uc002znh.2
chr22	ucsc	exon	18379490	18379747	.	-	.	ID=uc002znh.2:E13;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E13
chr22	ucsc	CDS	18379490	18379747	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS13;protein_id=uc002znh.2
chr22	ucsc	exon	18382214	18382314	.	-	.	ID=uc002znh.2:E14;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E14
chr22	ucsc	CDS	18382214	18382314	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS14;protein_id=uc002znh.2
chr22	ucsc	exon	18383608	18383763	.	-	.	ID=uc002znh.2:E15;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E15
chr22	ucsc	CDS	18383608	18383763	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS15;protein_id=uc002znh.2
chr22	ucsc	exon	18384644	18384745	.	-	.	ID=uc002znh.2:E16;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E16
chr22	ucsc	CDS	18384644	18384745	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS16;protein_id=uc002znh.2
chr22	ucsc	exon	18385397	18385513	.	-	.	ID=uc002znh.2:E17;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E17
chr22	ucsc	CDS	18385397	18385513	.	-	2	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS17;protein_id=uc002znh.2
chr22	ucsc	exon	18387398	18387605	.	-	.	ID=uc002znh.2:E18;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E18
chr22	ucsc	CDS	18387398	18387605	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS18;protein_id=uc002znh.2
chr22	ucsc	exon	18389315	18389652	.	-	.	ID=uc002znh.2:E19;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E19
chr22	ucsc	CDS	18389315	18389578	.	-	0	Parent=transcript:uc002znh.2.g1;ID=CDS:uc002znh.2:CDS19;protein_id=uc002znh.2
chr22	ucsc	exon	18507047	18507325	.	-	.	ID=uc002znh.2:E20;Parent=transcript:uc002znh.2.g1;Name=uc002znh.2;biotype=protein_coding;exon_id=uc002znh.2:E20
```

