# Biostar173114

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

make a bam file smaller by removing unwanted information see also https://www.biostars.org/p/173114/


## Usage

```
Usage: biostar173114 [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 9
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -keepAtt, --keepAttributes
      keep Attributes
      Default: false
    -keepCigar, --keepCigar
      keep cigar : don't remove hard clip
      Default: false
    -keepName, --keepName, --name
      keep Read Name, do not try to create a shorter name
      Default: false
    -keepQuals, --keepQualities
      keep base qualities
      Default: false
    -keepRG, --keepReadGroup
      if attributes are removed, keep the RG
      Default: false
    -keepSeq, --keepSequence
      keep read sequence
      Default: false
    -mate, --mate
      keep Mate/Paired Information
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
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
$ ./gradlew biostar173114
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar173114.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/biostar/Biostar173114.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar173114Test.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/biostar/Biostar173114Test.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **biostar173114** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



### Example

```
 $ java -jar dist/biostar173114.jar --keepSequence    my.bam  

@HD	VN:1.5	GO:none	SO:coordinate
@SQ	SN:rotavirus	LN:1074
R0	0	rotavirus	1	60	70M	*	0	0	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATT	*
R1	0	rotavirus	1	60	70M	*	0	0	GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATT	*
R2	0	rotavirus	1	60	70M	*	0	0	GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATT	*
R3	0	rotavirus	1	60	70M	*	0	0	GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATT	*
R4	0	rotavirus	1	60	70M	*	0	0	GGCATTTAATGCTTAACAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTCTTATT	*
R5	0	rotavirus	2	60	70M	*	0	0	GCTTTTAAAGCTTTTCAGTTGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGTTACGCTCTATTATTA	*
R6	0	rotavirus	2	60	51M19S	*	0	0	GCTTTTAATGCTTTTCAGTTGTTGCTGCACAAGATGGAGTCTACACAGCAGCTGTTCATCTCTCTTCATC	*
R7	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATGCTTTTCAGTGGTTTCTTCTCACGATGGAGTCTACTCAGCAGAAGGTAAGCACTATTATTA	*
R8	0	rotavirus	2	60	70M	*	0	0	GCTTTTAAAGCATTACAGTTGTTGCAGCTCAAGAAGGAGACTACTCAGCAGATGGTAAGCTCTATAATTA	*
R9	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATTCTATTCAGTGGTTGCTGCTCCAGAAGGAGTCTACTCAGGAGATGGTACGCTCTCTTATTA	*
Ra	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCAATTATTA	*
Rb	0	rotavirus	2	60	70M	*	0	0	GCTTTTAATGCTTTTCAGTTGTAGCTGCTCAAGATGGAGTCTACTCATCAGATGGTAAGCTCTCTTCTTA	*
Rc	0	rotavirus	2	60	63M7S	*	0	0	TCTTTAAATGCTTTTCAGTGTTTGCTGCTCAAGATGGAGTCTACTCAGCAGAAGGTAAGCTCTCTTAAAC	*
Rd	0	rotavirus	2	60	66M4S	*	0	0	GCATTTAATGCTTTTCAGTGGTTGCTGCACAAGATGGAGTCTACTCAGCAGATTGTAAGCTCTATTCTAA	*
Re	0	rotavirus	3	60	70M	*	0	0	CTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGAAGGCGTCTCCTGATGAGATGGTAAGCTCTATTATTAA	*
Rf	0	rotavirus	3	60	70M	*	0	0	CTTTTAATGGTTATGAGTGGTTGGTGCACAAGATGGAGTCTACTCAGCAGATGGTACTCTCTATAATTAA	*
R10	0	rotavirus	3	60	70M	*	0	0	CTTTTAAAGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTACTCTCTATTCTTAA	*
R11	0	rotavirus	3	60	70M	*	0	0	CTTTTAAAGCTTTTCAGAGGTTGCTGCTCAAGATGTAGTCTACTCAGGAGATTGTAAGCTCTATTATTAA	
```

## See also

  * https://bioinformatics.stackexchange.com/a/3866/71


