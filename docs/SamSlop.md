# SamSlop

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

extends sam by 'x' bases


## Usage

```
Usage: SamSlop [options] Files
  Options:
    --bamcompression
      Compression Level.
      Default: 5
    -q, --defaultQual
      default base quality (one letter)
      Default: #
    -M, --extend3
      num bases to extends on 3'
      Default: 0
    -m, --extend5
      num bases to extends on 5'
      Default: 0
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -r, --reference
      Indexed fasta Reference
    -c, --rmClip
      remove clipped bases
      Default: false
    --samoutputformat
      Sam output format.
      Default: SAM
      Possible Values: [BAM, SAM, CRAM]
    --version
      print version and exit

```

## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew SamSlop
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamSlop.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamSlop.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **SamSlop** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)





### Example



```

$ java -jar dist/samslop.jar -m 20 -M 10 -r ~/src/gatk-ui/testdata/ref.fa  ~/src/gatk-ui/testdata/S1.bam
 
@HD     VN:1.5  GO:none SO:unsorted
@SQ     SN:rotavirus    LN:1074
@RG     ID:S1   SM:S1
@PG     ID:bwa  PN:bwa  VN:0.7.12-r1044 CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_01_R1.fq.gz S1_01_R2.fq.gz
@PG     ID:bwa.1        PN:bwa  VN:0.7.12-r1044 CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_02_R1.fq.gz S1_02_R2.fq.gz
@PG     ID:bwa.2        PN:bwa  VN:0.7.12-r1044 CL:../bwa/bwa mem -R @RG\tID:S1\tSM:S1 ref.fa S1_03_R1.fq.gz S1_03_R2.fq.gz
rotavirus_1_317_5:0:0_7:0:0_2de 99      rotavirus       1       60      140M    =       248     317     GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAATATGGCGTCAACTCAGCAGATGGTCAGCTCTAATATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:33G4A3T14A7T4      RG:Z:S1 NM:i:5  AS:i:45 XS:i:0
rotavirus_1_535_4:0:0_4:0:0_1a6 163     rotavirus       1       60      140M    =       466     535     GGCTTTTACTGCTTTTCAGTGGTTGCTTCTCAAGATGGAGTGTACTCATCAGATGGTAAGCTCTATTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:8A18G13C6G21       RG:Z:S1 NM:i:4  AS:i:50 XS:i:0
rotavirus_1_543_5:0:0_11:0:0_390        163     rotavirus       1       60      140M    =       487     530     GGCTTTTAATGCTTTTCATTTGATGCTGCTCAAGATGGAGTCTACACAGCAGATGGTCAGCTCTATTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:18G1G1T22T11A12    RG:Z:S1 NM:i:5  AS:i:45      XS:i:0
rotavirus_1_578_3:0:0_7:0:0_7c  99      rotavirus       1       60      140M    =       509     578     GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTCCTGAGCAGCTGGTAAGCTCTATTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:43A2C5A17  RG:Z:S1 NM:i:3  AS:i:55 XS:i:0
rotavirus_1_497_4:0:0_5:0:0_2f6 163     rotavirus       1       60      140M    =       432     497     GGCATTTAATGCTTAACAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCTCTTATTAATACTTCTT     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########        MD:Z:3T10T0T48A5        RG:Z:S1 NM:i:4  AS:i:51 XS:i:0
rotavirus_2_509_4:0:0_8:0:0_68  99      rotavirus       1       60      4M      =       440     478     GGCTTTTAAAGCTTTTCAGTTGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGTTACGCTCTATTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:8T10G34G2A12       RG:Z:S1 NM:i:4  AS:i:50 XS:i:0
rotavirus_2_441_11:0:0_4:0:0_298        99      rotavirus       1       60      4M      =       372     440     GGCTTTTAATGCTTTTCAGTTGTTGCTGCACAAGATGGAGTCTACACAGCAGATGGTAAGCT       #+++++++++++++++++++++++++++++++++++++++++++++++++++##########  MD:Z:19G8T15T6  RG:Z:S1 NM:i:3  AS:i:36 XS:i:0
rotavirus_2_538_5:0:0_6:0:0_6b  99      rotavirus       1       60      4M      =       469     537     GGCTTTTAATGCTTTTCAGTGGTTTCTTCTCACGATGGAGTCTACTCAGCAGAAGGTAAGCACTATTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:23G2G4A20T7T9      RG:Z:S1 NM:i:5  AS:i:45 XS:i:0
rotavirus_2_459_8:0:0_3:0:0_388 99      rotavirus       1       60      4M      =       390     458     GGCTTTTAAAGCATTACAGTTGTTGCAGCTCAAGAAGGAGACTACTCAGCAGATGGTAAGCTCTATAATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:8T2T2T4G5T8T4T25T4 RG:Z:S1 NM:i:8  AS:i:30 XS:i:0
rotavirus_2_577_7:0:0_6:0:0_64  163     rotavirus       1       60      4M      =       513     576     GGCTTTTAATTCTATTCAGTGGTTGCTGCTCCAGAAGGAGTCTACTCAGGAGATGGTACGCTCTCTTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:9G2T17A3T13C8A5A6  RG:Z:S1 NM:i:7  AS:i:35 XS:i:0
rotavirus_2_500_1:0:0_14:0:0_88 163     rotavirus       1       60      4M      =       453     492     GGCTTTTAATGCTTTTCAGTGGTTGCTGCTCAAGATGGAGTCTACTCAGCAGATGGTAAGCTCAATTATTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:62T7       RG:Z:S1 NM:i:1  AS:i:65 XS:i:0
rotavirus_2_445_5:0:0_10:0:0_c9 163     rotavirus       1       60      4M      =       377     434     GGCTTTTAATGCTTTTCAGTTGTAGCTGCTCAAGATGGAGTCTACTCATCAGATGGTAAGCTCTCTTCTTAATACTTCTTT    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########       MD:Z:19G2T24G15A2A3     RG:Z:S1 NM:i:5  AS:i:48 XS:i:0
rotavirus_2_488_8:0:0_9:0:0_343 99      rotavirus       1       60      4M      =       419     487     GTCTTTAAATGCTTTTCAGTGTTTGCTGCTCAAGATGGAGTCTACTCAGCAGAAGGTAAGCTCTATTATTAATA   #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########      MD:Z:0G4T14G31T10       RG:Z:S1 NM:i:4  AS:i:47 XS:i:0
rotavirus_2_478_5:0:0_5:0:0_35c 163     rotavirus       1       60      4M      =       409     477     GGCATTTAATGCTTTTCAGTGGTTGCTGCACAAGATGGAGTCTACTCAGCAGATTGTAAGCTCTATTATTAATACTT        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##########   MD:Z:2T25T24G12 RG:Z:S1 NM:i:3  AS:i:53 XS:i:0
(...)
```



