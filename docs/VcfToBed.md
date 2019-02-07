# VcfToBed

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

vcf to bed


## Usage

```
Usage: vcf2bed [options] Files
  Options:
    -F, --format
      output format
      Default: bed
      Possible Values: [bed, interval]
    -header, --header
      Print Header
      Default: false
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -M, --max
      Optional filter: max sequence length. A distance specified as a positive 
      integer.Comma are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
    -m, --min
      Optional filter: min sequence length. A distance specified as a positive 
      integer.Comma are removed. The following suffixes are interpreted : 
      b,bp,k,kb,m,mb 
    -c, --no-ci
      For structural variant, ignore the extention of the boundaries using 
      INFO/CIPOS and INFO/CIEND
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -R, --reference
      Convert the contigs of the VCF on the fly using an indexed genome. The 
      parameter is the path to an Indexed fasta Reference file. This fasta 
      file must be indexed with samtools faidx and with picard 
      CreateSequenceDictionary. The parameter can also be a 'key' (matching 
      the regular expression `[A-Za-z][A-Za-z0-9_\\-]*`) in a catalog file. A 
      'catalog' file is a java property file ( 
      https://docs.oracle.com/javase/tutorial/essential/environment/properties.html 
      ) where the values are the path to the fasta file.  Catalogs are 
      searched in that order : `${PWD}/fasta-ref.properties`, 
      `${HOME}/.fasta-ref.properties`, `/etc/jvarkit/fasta-ref.properties`.  
      If the key or the path are not defined by the user, they will be 
      searched in that order 1) the java property 
      -Djvarkit.fasta.reference=pathTofastaOrCatalogKey . 2) the linux 
      environement variable $FASTA_REFERENCE=pathTofastaOrCatalogKey 3) The 
      catalogs. 
    -x, --slop
      Extends interval by 'x' bases on both sides. A distance specified as a 
      positive integer.Comma are removed. The following suffixes are 
      interpreted : b,bp,k,kb,m,mb
      Default: 0
    --version
      print version and exit

```


## Keywords

 * bed
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcf2bed
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToBed.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VcfToBed.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfToBedTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VcfToBedTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcf2bed** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Motivation

I'm lazy about using awk or bioalcidaejdk for this task and I want something that uses INFO/CIPOS and INFO/CIEND for structural variants

## Input

input is one or more VCF file

one file ending with '.list' is interpreted as a list of paths (one per lines)

if there is no input, the program reads vcf from stdin

##Example


```
$ wget -q -O - "https://github.com/hall-lab/cshl_sv_2014/blob/master/supplemental/NA12878.lumpy.vcf?raw=true" |\
	grep -A 10 '#CHROM'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	869423	1	G	<DEL>	345.50	.	SVTYPE=DEL;SVLEN=-857;END=870280;STR=+-:25;IMPRECISE;CIPOS=-1,34;CIEND=0,0;EVENT=1;SUP=25;PESUP=25;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	1/1:25:25:0:0.00:71:19:52:345.50:-36,-8,-1
1	1588585	5	A	<DUP>	0.00	.	SVTYPE=DUP;SVLEN=65356;END=1653941;STR=-+:7;IMPRECISE;CIPOS=-126,1;CIEND=-2,67;EVENT=5;SUP=7;PESUP=7;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:7:7:0:74.40:139:125:13:0.00:-1,-8,-24
1	1594964	6	C	<DUP>	0.00	.	SVTYPE=DUP;SVLEN=65855;END=1660819;STR=-+:8;IMPRECISE;CIPOS=-81,2;CIEND=-1,127;EVENT=6;SUP=8;PESUP=8;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:8:8:0:77.96:153:137:15:0.00:-1,-9,-25
1	2566176	7	A	<DEL>	121.20	.	SVTYPE=DEL;SVLEN=-418;END=2566594;STR=+-:14;IMPRECISE;CIPOS=-2,68;CIEND=0,0;EVENT=7;SUP=14;PESUP=14;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:14:14:0:0.00:78:44:33:121.20:-13,-1,-12
1	2911548	8	G	<DEL>	440.34	.	SVTYPE=DEL;SVLEN=-302;END=2911850;STR=+-:20;CIPOS=0,0;CIEND=0,0;EVENT=8;SUP=20;PESUP=8;SRSUP=12;EVTYPE=PE,SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:20:8:12:0.00:185:86:99:440.34:-48,-4,-15
1	2919034	9	G	<DEL>	289.83	.	SVTYPE=DEL;SVLEN=-332;END=2919366;STR=+-:22;CIPOS=0,0;CIEND=0,0;EVENT=9;SUP=22;PESUP=10;SRSUP=12;EVTYPE=PE,SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:22:10:12:0.00:160:86:74:289.83:-31,-2,-20
1	5447229	14	G	<DUP>	380.12	.	SVTYPE=DUP;SVLEN=210;END=5447439;STR=-+:11;CIPOS=0,0;CIEND=0,0;EVENT=14;SUP=11;PESUP=1;SRSUP=10;EVTYPE=PE,SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	1/1:11:1:10:0.00:197:104:93:380.12:-39,-7,-1
1	5876603	15	G	<DEL>	0.00	.	SVTYPE=DEL;SVLEN=-928;END=5877531;STR=+-:8;CIPOS=0,0;CIEND=0,0;EVENT=15;SUP=8;PESUP=0;SRSUP=8;EVTYPE=SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:8:0:8:294.07:169:168:0:0.00:-8,-37,-117
1	5877530	16	T	<DEL>	63.31	.	SVTYPE=DEL;SVLEN=-72;END=5877602;STR=+-:13;CIPOS=0,0;CIEND=0,0;EVENT=16;SUP=13;PESUP=0;SRSUP=13;EVTYPE=SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:13:0:13:0.00:188:136:51:63.31:-10,-4,-54
1	6619067	19_1	T	[1:6619506[T	0.00	.	SVTYPE=BND;STR=--:7;IMPRECISE;CIPOS=-88,1;CIEND=-26,2;MATEID=19_2;EVENT=19;SUP=7;PESUP=7;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:7:7:0:127.76:131:117:13:0.00:-1,-14,-66


$ wget -q -O - "https://github.com/hall-lab/cshl_sv_2014/blob/master/supplemental/NA12878.lumpy.vcf?raw=true" |\
	java -jar dist/vcf2bed.jar |\
	head

1	869421	870280	1	345
1	1588458	1654008	5	0
1	1594882	1660946	6	0
1	2566173	2566594	7	121
1	2911547	2911850	8	440
1	2919033	2919366	9	289
1	5447228	5447439	14	380
1	5876602	5877531	15	0
1	5877529	5877602	16	63
1	6618978	6619069	19_1	0
```

