# HowManyBamDict

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

finds if there's are some differences in the sequence dictionaries.


## Usage

```
Usage: howmanybamdict [options] Files
  Options:
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

 * sam
 * bam
 * dict


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew howmanybamdict
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/HowManyBamDict.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/HowManyBamDict.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **howmanybamdict** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ find PROJ/ -name "accepted_hits.bam" |\
 java -jar dist/howmanybamdict.jar

DICT	208e057b92b5c45262c59b40209b5fde	22	2725537669	chr1=195471971;chr10=130694993;chr11=122082543;chr12=120129022;chr13=120421639;chr14=124902244;chr15=104043685;chr16=98207768;chr17=94987271;chr18=90702639;chr19=61431566;chr2=182113224;chr3=160039680;chr4=156508116;chr5=151834684;chr6=149736546;chr7=145441459;chr8=129401213;chr9=124595110;chrM=16299;chrX=171031299;chrY=91744698	PROJ/S1/accepted_hits.bams.bam
BAM	PROJ/S1/accepted_hits.bam	208e057b92b5c45262c59b40209b5fde
DICT	86b39aef44f740da797b89baf3f505d8	66	2730871774	chr1=195471971;chr10=130694993;chr11=122082543;chr12=120129022;chr13=120421639;chr14=124902244;chr15=104043685;chr16=98207768;chr17=94987271;chr18=90702639;chr19=61431566;chr1_GL456210_random=169725;chr1_GL456211_random=241735;chr1_GL456212_random=153618;chr1_GL456213_random=39340;chr1_GL456221_random=206961;chr2=182113224;chr3=160039680;chr4=156508116;chr4_GL456216_random=66673;chr4_GL456350_random=227966;chr4_JH584292_random=14945;chr4_JH584293_random=207968;chr4_JH584294_random=191905;chr4_JH584295_random=1976;chr5=151834684;chr5_GL456354_random=195993;chr5_JH584296_random=199368;chr5_JH584297_random=205776;chr5_JH584298_random=184189;chr5_JH584299_random=953012;chr6=149736546;chr7=145441459;chr7_GL456219_random=175968;chr8=129401213;chr9=124595110;chrM=16299;chrUn_GL456239=40056;chrUn_GL456359=22974;chrUn_GL456360=31704;chrUn_GL456366=47073;chrUn_GL456367=42057;chrUn_GL456368=20208;chrUn_GL456370=26764;chrUn_GL456372=28664;chrUn_GL456378=31602;chrUn_GL456379=72385;chrUn_GL456381=25871;chrUn_GL456382=23158;chrUn_GL456383=38659;chrUn_GL456385=35240;chrUn_GL456387=24685;chrUn_GL456389=28772;chrUn_GL456390=24668;chrUn_GL456392=23629;chrUn_GL456393=55711;chrUn_GL456394=24323;chrUn_GL456396=21240;chrUn_JH584304=114452;chrX=171031299;chrX_GL456233_random=336933;chrY=91744698;chrY_JH584300_random=182347;chrY_JH584301_random=259875;chrY_JH584302_random=155838;chrY_JH584303_random=158099	PROJ/S2/accepted_hits.bam	ROJ/S2/accepted_hits.bam
BAM	PROJ/S2/accepted_hits.bam	86b39aef44f740da797b89baf3f505d8
BAM	PROJ/S3/accepted_hits.bam	86b39aef44f740da797b89baf3f505d8
```

