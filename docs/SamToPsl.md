# SamToPsl

Convert SAM/BAM to PSL http://genome.ucsc.edu/FAQ/FAQformat.html#format2 or BED12


## Usage

```
Usage: sam2psl [options] Files
  Options:
    -h, --help
      print help and exit
    --helpFormat
      What kind of help
      Possible Values: [usage, markdown, xml]
    -o, --output
      Output file. Optional . Default: stdout
    -s, --single
      treat all reads as single end
      Default: false
    --version
      print version and exit
    -B, bed12
      export as BED 12
      Default: false

```

## Compilation

### Requirements / Dependencies

* java compiler SDK 1.8 http://www.oracle.com/technetwork/java/index.html (**NOT the old java 1.7 or 1.6**) . Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )
* GNU Make >= 3.81
* curl/wget
* git
* xsltproc http://xmlsoft.org/XSLT/xsltproc2.html (tested with "libxml 20706, libxslt 10126 and libexslt 815")


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ make sam2psl
```

The *.jar libraries are not included in the main jar file, so you shouldn't move them (https://github.com/lindenb/jvarkit/issues/15#issuecomment-140099011 ).
The required libraries will be downloaded and installed in the `dist` directory.

### edit 'local.mk' (optional)

The a file **local.mk** can be created edited to override/add some definitions.

For example it can be used to set the HTTP proxy:

```
http.proxy.host=your.host.com
http.proxy.port=124567
```
## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamToPsl.java
](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/SamToPsl.java
)
## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **sam2psl** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)



Convert **SAM/BAM** to **PSL** http://genome.ucsc.edu/FAQ/FAQformat.html#format2 or **BED12** .

Properly-paired reads are extended to the mate's position.

What ? **bamtobed** http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html does the same job ?! too late.




### Example




```
$ samtools view -b  http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeq10t12C3hFR2x75Th131Il200AlnRep1.bam "chr15:81575506-81616397" |\
   java -jar ~/src/jvarkit-git/dist/sam2psl.jar -s   > out.psl

$ tail out.psl
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1105:11756:193850_2:N:0:/2_147	101	0	101	chr15	10349497481616327	81622456	3	1,5,95,	0,1,6,	81616327,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1301:4643:94800_2:N:0:/2_147	101	0	101	chr15	103494974	81616334	81622456	3	1,5,95,	0,1,6,	81616334,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1308:18580:185579_2:Y:0:/2_147	101	0	101	chr15	10349497481616327	81622456	3	1,5,95,	0,1,6,	81616327,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2205:10117:76559_2:N:0:/2_147	101	0	101	chr15	10349497481616321	81622456	3	1,5,95,	0,1,6,	81616321,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2206:7885:15613_2:N:0:/2_403	101	0	101	chr15	103494974	81613633	81622456	3	1,5,95,	0,1,6,	81613633,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2206:7885:15613_2:N:0:/2_147	101	0	101	chr15	103494974	81616308	81622456	3	1,5,95,	0,1,6,	81616308,81616393,81622362,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:2303:12879:149117_1:Y:0:/1_83	101	0	101	chr15	10349497481616334	81622456	3	1,5,95,	0,1,6,	81616334,81616393,81622362,
6064	0	0	0	0	0	0	5964	+	HWI-ST0787:100:C02F9ACXX:3:1301:11600:100190_1:N:0:/1_99	100	0	100	chr15	10349497481616393	81622457	2	4,96,	0,4,	81616393,81622361,
6064	0	0	0	0	0	0	5964	+	HWI-ST0787:100:C02F9ACXX:3:2304:5980:187674_1:Y:0:/1_99	100	0	100	chr15	103494974	81616393	81622457	2	4,96,	0,4,	81616393,81622361,
6065	0	0	0	0	0	0	5964	-	HWI-ST0787:100:C02F9ACXX:3:1306:18607:99733_2:N:0:/2_147	101	0	101	chr15	10349497481616334	81622457	3	1,4,96,	0,1,5,	81616334,81616394,81622362,


```



used as a custom track in the **UCSC genome browser**.

![img](http://i.imgur.com/Gi6Sd0M.png)




### See also


* bedtools/bamtobed : http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html






