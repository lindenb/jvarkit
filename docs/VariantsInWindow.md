# VariantsInWindow

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Annotate Number of Variants overlaping a sliding window.


## Usage

```
Usage: variantsinwindow [options] Files
  Options:
    --best, -best
      Only print the window with the hightest number of matches
      Default: false
    -filter, --filter
      if --treshold is != -1 and the number of matches is greater than 
      threshold, set this FILTER
      Default: TOO_MANY_CLOSE_VARIANTS
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -o, --output
      Output file. Optional . Default: stdout
    -treshold, --treshold
      Number of variants to set the FILTER
      Default: -1
    -vf, --variant-filter
      Variants we want to keep. Variant FAILING that Jexl expression will be 
      excluded from the window.A Java EXpression Language (JEXL) expressions 
      to filter the variants from a VCF. Empty string will accept all 
      variants. Expression returning a TRUE will accept the variant. See 
      https://gatkforums.broadinstitute.org/gatk/discussion/1255 
      Default: <empty string> (ACCEPT ALL)
    --version
      print version and exit
    -S, --shift, --windowShift
      Window shift.A distance specified as a positive integer.Comma are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 50
    -W, --windowSize
      Window Size.A distance specified as a positive integer.Comma are 
      removed. The following suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 150
    -noemptywin
      Don't print Windows in INFO having zero match.
      Default: false
    -windowName
      INFO Attribute name that will be added
      Default: WINDOW

```


## Keywords

 * vcf
 * annotation



## See also in Biostars

 * [https://www.biostars.org/p/291144](https://www.biostars.org/p/291144)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew variantsinwindow
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VariantsInWindow.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/misc/VariantsInWindow.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VariantsInWindowTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/misc/VariantsInWindowTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **variantsinwindow** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Example

```
##fileformat=VCFv4.2
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=TOO_MANY_CLOSE_VARIANTS,Description="Filter defined in vcfwindowvariants">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
[lindenb@kaamelot-master01 jvarkit-git]$ java -jar dist/variantsinwindow.jar ~/src/gatk-ui/testdata/mutations.vcf --treshold 1 -shift 1 -windowSize 10
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SGB=-61.9012;VDB=3.53678e-05;WINDOW=41|50|1|0,42|51|1|0,43|52|1|0,44|53|1|0,45|54|1|0,46|55|1|0,47|56|1|0,48|57|1|0,49|58|1|0,50|59|1|0,51|60|1|0	GT:PL	0/0:0,255,134	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SGB=-369.163;VDB=0.313337;WINDOW=81|90|1|0,82|91|1|0,83|92|1|0,84|93|1|0,85|94|1|0,86|95|1|0,87|96|1|0,88|97|1|0,89|98|1|0,90|99|1|0,91|100|1|0	GT:PL	0/0:0,255,133	0/1:40,0,31	0/0:0,255,134	0/0:0,255,82
rotavirus	130	.	T	C	4.12	.	AC1=1;AF1=0.124933;BQB=1;DP=1349;DP4=1139,0,204,0;FQ=4.48321;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.762964;SGB=-335.275;VDB=0.00084636;WINDOW=120|129|1|0,121|130|1|0,122|131|1|0,123|132|1|0,124|133|1|0,125|134|1|0,126|135|1|0,127|136|1|0,128|137|1|0,129|138|1|0,130|139|1|0	GT:PL	0/1:38,0,35	0/0:0,255,132	0/0:0,255,132	0/0:0,255,79
rotavirus	232	.	T	A	5.45	.	AC1=1;AF1=0.124959;BQB=1;DP=1308;DP4=1098,0,207,0;FQ=5.87117;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.453119;SGB=-340.116;VDB=0.0100544;WINDOW=222|231|1|0,223|232|1|0,224|233|1|0,225|234|1|0,226|235|1|0,227|236|1|0,228|237|1|0,229|238|1|0,230|239|1|0,231|240|1|0,232|241|1|0	GT:PL	0/1:40,0,35	0/0:0,255,135	0/0:0,255,132	0/0:0,255,81
rotavirus	267	.	C	G	4.76	.	AC1=1;AF1=0.124951;BQB=0.953186;DP=1393;DP4=1156,0,234,0;FQ=5.15288;MQ=60;MQ0F=0;MQB=1;PV4=1,5.65123e-05,1,1;RPB=0.0284076;SGB=-383.686;VDB=0.367507;WINDOW=257|266|1|0,258|267|1|0,259|268|1|0,260|269|1|0,261|270|1|0,262|271|1|0,263|272|1|0,264|273|1|0,265|274|1|0,266|275|1|0,267|276|1|0	GT:PL	0/1:39,0,32	0/0:0,255,132	0/0:0,255,136	0/0:0,255,76
rotavirus	424	.	A	G	52.99	.	AC1=1;AF1=0.125;BQB=0.956333;DP=1555;DP4=1096,206,198,55;FQ=53.5713;MQ=60;MQ0F=0;MQB=1;MQSB=1;PV4=0.0270045,4.0796e-05,1,1;RPB=0.0600948;SGB=-160.094;VDB=0.000623759;WINDOW=414|423|1|0,415|424|1|0,416|425|1|0,417|426|1|0,418|427|1|0,419|428|1|0,420|429|1|0,421|430|1|0,422|431|1|0,423|432|1|0,424|433|1|0	GT:PL	0/0:0,255,200	0/1:89,0,51	0/0:0,255,189	0/0:0,255,153
rotavirus	520	.	T	A	53.99	.	AC1=1;AF1=0.125;BQB=0.215002;DP=2372;DP4=1055,856,223,231;FQ=54.5713;MQ=60;MQ0F=0;MQB=1;MQSB=1;PV4=0.0211343,5.03431e-23,1,1;RPB=0.805995;SGB=-738.702;VDB=0.730588;WINDOW=510|519|1|0,511|520|1|0,512|521|1|0,513|522|1|0,514|523|1|0,515|524|1|0,516|525|1|0,517|526|1|0,518|527|1|0,519|528|1|0,520|529|1|0GT:PL	0/0:0,255,231	0/0:0,255,225	0/1:90,0,112	0/0:0,255,204
(...)
rotavirus	1054	.	C	G	15.65	TOO_MANY_CLOSE_VARIANTS	AC1=2;AF1=0.249999;BQB=1;DP=487;DP4=0,364,0,120;FQ=16.8692;G3=0.75,2.21169e-28,0.25;HWE=0.0339211;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.95941;SGB=42.7815;VDB=1.4013e-45;WINDOW=1044|1053|3|0,1045|1054|2|0,1046|1055|1|0,1047|1056|1|0,1048|1057|1|0,1049|1058|1|0,1050|1059|1|0,1051|1060|1|0,1052|1061|1|0,1053|1062|1|0,1054|1063|2|0	GT:PL	0/0:0,255,90	1/1:63,235,0	0/0:0,255,99	0/0:0,132,66
rotavirus	1064	.	G	A	21.56	TOO_MANY_CLOSE_VARIANTS	AC1=2;AF1=0.25;BQB=0.683886;DP=250;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16;WINDOW=1054|1063|2|0,1055|1064|1|0,1056|1065|1|0,1057|1066|1|0,1058|1067|1|0,1059|1068|1|0,1060|1069|1|0,1061|1070|1|0,1062|1071|1|0,1063|1072|1|0,1064|1073|1|0	GT:PL	0/0:0,244,70	0/0:0,199,65	0/0:0,217,68	1/1:69,84,0
```

