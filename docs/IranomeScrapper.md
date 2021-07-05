# IranomeScrapper

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Iranome scrapper


## Usage

```
Usage: iranomescrapper [options] Files
  Options:
    -failed, --failed
      Save the variants where no data in iranome data was found in this file.
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -T, --tabix
      Existing tabix-indexed database. Don't lookup the data if the variant is 
      already in this database.
    --version
      print version and exit

```


## Keywords

 * iranome
 * vcf


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew iranomescrapper
```

The java jar file will be installed in the `dist` directory.


## Creation Date

20210728

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/iranome/IranomeScrapper.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/iranome/IranomeScrapper.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **iranomescrapper** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
$ java -jar dist/iranomescrapper.jar src/test/resources/test_vcf01.vcf
[WARN][IranomeScrapper]Too many requests for http://www.iranome.ir/variant/1-912049-T-C max_sleep_millisec=1
[WARN][IranomeScrapper]Too many requests for http://www.iranome.ir/variant/1-913889-G-A max_sleep_millisec=1001
CHROM	POS	REF	ALT	AC	AN	AF	AC_ARAB	AC_AZERI	AC_BALOCH	AC_KURD	AC_LUR	AC_PERSIAN	AC_PERSIAN_GULF_ISLANDER	AC_TURKMEN	AN_ARAB	AN_AZERI	AN_BALOCH	AN_KURD	AN_LUR	AN_PERSIAN	AN_PERSIAN_GULF_ISLANDER	AN_TURKMEN	date
1	909238	G	C	1241	1586	0.782472	160	149	171	151	156	158	157	139	200	198	200	200	200	200	200	188	20210629_103622
1	912049	T	C	856	1558	0.549422	104	106	120	85	105	117	116	103	190	196	194	194	196	200	194	194	20210629_103623
1	914333	C	G	961	1584	0.606692	119	117	132	107	121	127	129	109	200	200	200	200	200	200	200	184	20210629_103626
1	914852	G	C	939	1596	0.588346	117	116	131	99	112	122	127	115	200	200	200	200	200	200	200	196	20210629_103628
1	914940	T	C	937	1598	0.586358	116	116	131	99	112	122	127	114	200	200	200	200	200	200	200	198	20210629_103628
1	935222	C	A	989	1554	0.636422	123	121	139	108	127	130	139	102	200	198	200	200	200	200	200	156	20210629_103638
1	949608	G	A	460	1600	0.2875	55	59	60	62	50	59	56	59	200	200	200	200	200	200	200	200	20210629_103647
1	984302	T	C	934	1554	0.60103	116	120	122	115	125	124	125	87	200	200	200	198	200	200	200	156	20210629_103719
1	985266	C	T	934	1600	0.58375	115	118	119	108	126	117	124	107	200	200	200	200	200	200	200	200	20210629_103720
1	985446	G	T	647	1330	0.486466	79	91	97	70	81	87	99	43	166	182	186	176	168	158	174	120	20210629_103721
1	985450	G	A	257	1352	0.190089	33	36	42	24	30	24	38	30	166	176	180	176	172	172	176	134	20210629_103723
1	1007432	G	A	908	1598	0.56821	114	111	109	115	113	115	123	108	200	200	200	200	200	200	200	198	20210629_103735
1	1018144	T	C	907	1596	0.568296	115	111	108	115	113	115	127	103	200	200	200	200	200	200	200	196	20210629_103746
1	1019180	T	C	947	1600	0.591875	116	114	113	118	118	119	131	118	200	200	200	200	200	200	200	200	20210629_103750
```

