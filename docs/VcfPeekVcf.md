# VcfPeekVcf

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Get the INFO from a VCF and use it for another VCF


## Usage

```
Usage: vcfpeekvcf [options] Files
  Options:
    -a, -alt, --alt
      How alt allele must be found in the variants of the indexed file.
      Default: none
      Possible Values: [none, all, at_least_one]
    -b, --buffer-size
      buffer size (in bp). We don't do a random access for each variant. 
      Instead of this, load all the variants in a defined window. A distance 
      specified as a positive integer.Comma are removed. The following 
      suffixes are interpreted : b,bp,k,kb,m,mb
      Default: 100000
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    -missingIsError, --missingIsError
      Missing Info Header is an error
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -p, --prefix
      prefix all database tags with this prefix to avoid collisions
      Default: <empty string>
    -i, --replaceid
      Replace the ID field if it exists
      Default: false
  * -f, --tabix, --resource
      The VCF file indexed with TABIX or tribble. Source of the annotations
    -span, --span
      [20180713] when checking for the '--alt' option, ignore spanning 
      deletion: *
      Default: false
    -t, --tags
      tag1,tag2,tag... the INFO keys to peek from the indexed file
      Default: []
    --version
      print version and exit

```


## Keywords

 * vcf
 * annotation


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcfpeekvcf
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcf.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcf.java)

### Unit Tests

[https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcfTest.java](https://github.com/lindenb/jvarkit/tree/master/src/test/java/com/github/lindenb/jvarkit/tools/vcfvcf/VcfPeekVcfTest.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcfpeekvcf** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

[http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)

> Lindenbaum, Pierre (2015): JVarkit: java-based utilities for Bioinformatics. figshare.
> [http://dx.doi.org/10.6084/m9.figshare.1425030](http://dx.doi.org/10.6084/m9.figshare.1425030)


## Alternate tools

you can also use `GATK VariantAnnotator` or `bcftools`. But this tool contains some interesting options.

## Example

Peek GENEINFO and VP from dbsnp_135.b37 to decorate /ExAC.r0.3.sites.vep.vcf.gz


```bash
$ java -jar dist/vcfpeekvcf.jar -f dbsnp_135.b37.vcf.gz -t GENEINFO,VP -p 00_NCBI135_  -i ExAC.r0.3.sites.vep.vcf.gz |\
grep NCBI135_ 


##INFO=<ID=00_NCBI135_GENEINFO,Number=1,Type=String,Description="Pairs each of gene symbol:gene id.  The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)">
##INFO=<ID=00_NCBI135_VP,Number=1,Type=String,Description="Variation Property">
##VcfPeekVcfCmdLine=-f /commun/data/pubdb/broadinstitute.org/bundle/1.5/b37/dbsnp_135.b37.vcf.gz -t GENEINFO,VP -p 00_NCBI135_ -i /commun/data/pubdb/broadinstitute.org/exac/0.3/ExAC.r0.3.sites.vep.vcf.gz
1	69428	rs140739101	T	G	4604953.42	PASS	00_NCBI135_GENEINFO=.;00_NCBI135_VP=050200000004040000000100;AC=2141;AC_AFR=20;AC_AMR=53;AC_Adj=1985;AC_EAS=0;AC_FIN=166;AC_Het=493;AC_Hom=746;AC_NFE=1668;AC_OTH=17;AC_SAS=61;AF=0.022;AN=99358;AN_AFR=7834;AN_AMR=6588;AN_Adj=80618;AN_EAS=8396;AN_FIN=3590;AN_NFE=41104;AN_OTH=626;AN_SAS=12480;BaseQRankSum=2.19;ClippingRankSum=-4.240e-01;DB;DP=1131603;FS=11.924;GQ_HIST=1523|6341|94|79|1292|32|35|13|11|10|20|12|24526|8836|2118|1909|782|234|274|1538,7|17|25|10|7|6|16|4|9|5|15|10|6|2|7|4|7|1|5|1190;GQ_MEAN=75.91;GQ_STDDEV=202.08;Het_AFR=10;Het_AMR=11;Het_EAS=0;Het_FIN=28;Het_NFE=418;Het_OTH=5;Het_SAS=21;Hom_AFR=5;Hom_AMR=21;Hom_EAS=0;Hom_FIN=69;Hom_NFE=625;Hom_OTH=6;Hom_SAS=20;InbreedingCoeff=0.3731;MQ=27.17;MQ0=0;MQRankSum=-1.014e+00;NCC=18661;QD=14.64;ReadPosRankSum=1.08;VQSLOD=44.51;culprit=MQ
1	69453	rs142004627	G	A	6724.79	VQSRTrancheSNP99.90to99.95	00_NCBI135_GENEINFO=.;00_NCBI135_VP=050200000004000000000100;AC=12;AC_AFR=3;AC_AMR=4;AC_Adj=7;AC_EAS=0;AC_FIN=0;AC_Het=1;AC_Hom=3;AC_NFE=0;AC_OTH=0;AC_SAS=0;AF=1.181e-04;AN=101644;AN_AFR=7954;AN_AMR=6638;AN_Adj=81032;AN_EAS=8402;AN_FIN=3606;AN_NFE=41318;AN_OTH=626;AN_SAS=12488;BaseQRankSum=0.742;ClippingRankSum=-6.150e-01;DB;DP=1001026;FS=90.555;GQ_HIST=1494|6818|211|263|1616|53|26|16|8|3|4|3|25002|9059|2178|1982|821|251|319|695,1|0|0|0|1|1|1|0|1|0|0|0|0|0|0|0|1|0|0|2;GQ_MEAN=53.93;GQ_STDDEV=24.34;Het_AFR=1;Het_AMR=0;Het_EAS=0;Het_FIN=0;Het_NFE=0;Het_OTH=0;Het_SAS=0;Hom_AFR=1;Hom_AMR=2;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;InbreedingCoeff=0.0192;MQ=28.57;MQ0=0;MQRankSum=1.54;NCC=17104;QD=22.95;ReadPosRankSum=-4.060e-01;VQSLOD=-4.645e+01;culprit=FS
1	69496	rs150690004	G	A	323905.92	VQSRTrancheSNP99.60to99.80	00_NCBI135_GENEINFO=.;00_NCBI135_VP=050200000004040000000100;AC=62;AC_AFR=46;AC_AMR=11;AC_Adj=59;AC_EAS=0;AC_FIN=0;AC_Het=35;AC_Hom=12;AC_NFE=2;AC_OTH=0;AC_SAS=0;AF=6.729e-04;AN=92132;AN_AFR=7748;AN_AMR=6566;AN_Adj=79928;AN_EAS=8388;AN_FIN=3570;AN_NFE=40582;AN_OTH=624;AN_SAS=12450;BaseQRankSum=2.12;ClippingRankSum=-1.529e+00;DB;DP=991114;FS=8.274;GQ_HIST=893|4697|202|135|905|245|364|326|147|307|231|164|22619|8691|2105|1925|794|250|316|750,0|0|0|0|1|0|0|0|0|1|1|0|0|0|0|1|0|0|0|46;GQ_MEAN=57.92;GQ_STDDEV=73.03;Het_AFR=24;Het_AMR=9;Het_EAS=0;Het_FIN=0;Het_NFE=2;Het_OTH=0;Het_SAS=0;Hom_AFR=11;Hom_AMR=1;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;InbreedingCoeff=0.0180;MQ=40.87;MQ0=0;MQRankSum=1.84;NCC=23714;NEGATIVE_TRAIN_SITE;QD=32.93;ReadPosRankSum=2.93;VQSLOD=-3.741e+00;culprit=MQ
1	69511	rs75062661	A	G	120729371.20	PASS	00_NCBI135_GENEINFO=OR4F5:79501;00_NCBI135_VP=050200000000000110000100;AC=75589;AC_AFR=4392;AC_AMR=6155;AC_Adj=72743;AC_EAS=8379;AC_FIN=3289;AC_Het=1789;AC_Hom=35477;AC_NFE=37731;AC_OTH=572;AC_SAS=12225;AF=0.894;AN=84570;AN_AFR=7392;AN_AMR=6474;AN_Adj=77432;AN_EAS=8384;AN_FIN=3320;AN_NFE=38832;AN_OTH=596;AN_SAS=12434;BaseQRankSum=0.831;ClippingRankSum=1.06;DB;DP=3157075;FS=23.500;GQ_HIST=1057|1461|500|275|204|94|127|139|117|163|175|204|1594|421|305|431|495|416|553|33554,71|556|484|263|161|89|124|138|114|159|175|198|273|346|295|424|491|416|553|33550;GQ_MEAN=224.54;GQ_STDDEV=255.92;Het_AFR=1030;Het_AMR=157;Het_EAS=5;Het_FIN=11;Het_NFE=461;Het_OTH=12;Het_SAS=113;Hom_AFR=1681;Hom_AMR=2999;Hom_EAS=4187;Hom_FIN=1639;Hom_NFE=18635;Hom_OTH=280;Hom_SAS=6056;InbreedingCoeff=0.6382;MQ=31.34;MQ0=0;MQRankSum=-4.020e-01;NCC=29303;QD=26.34;ReadPosRankSum=-1.106e+00;VQSLOD=131.28;culprit=FS
1	69534	rs190717287	T	C	109944.08	PASS	00_NCBI135_GENEINFO=.;00_NCBI135_VP=050200000000000010000100;AC=27;AC_AFR=0;AC_AMR=0;AC_Adj=26;AC_EAS=26;AC_FIN=0;AC_Het=20;AC_Hom=3;AC_NFE=0;AC_OTH=0;AC_SAS=0;AF=3.005e-04;AN=89844;AN_AFR=7810;AN_AMR=6532;AN_Adj=78908;AN_EAS=8392;AN_FIN=3344;AN_NFE=39788;AN_OTH=606;AN_SAS=12436;BaseQRankSum=5.94;ClippingRankSum=0.131;DB;DP=1383040;FS=0.000;GQ_HIST=424|4210|191|117|539|33|26|18|3|4|1|2|28777|7095|1331|1167|439|98|77|370,0|0|0|0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|23;GQ_MEAN=58.29;GQ_STDDEV=53.52;Het_AFR=0;Het_AMR=0;Het_EAS=20;Het_FIN=0;Het_NFE=0;Het_OTH=0;Het_SAS=0;Hom_AFR=0;Hom_AMR=0;Hom_EAS=3;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;InbreedingCoeff=0.0081;MQ=38.81;MQ0=0;MQRankSum=0.393;NCC=26044;QD=16.32;ReadPosRankSum=0.777;VQSLOD=-1.477e+00;culprit=MQ
1	69552	rs55874132	G	T,A,C	6289.25	VQSRTrancheSNP99.60to99.80	00_NCBI135_GENEINFO=OR4F5:79501;00_NCBI135_VP=050300000000040400000100;AC=3,3,5;AC_AFR=0,0,0;AC_AMR=3,0,0;AC_Adj=3,3,0;AC_EAS=0,0,0;AC_FIN=0,0,0;AC_Het=1,1,0,0,0,0;AC_Hom=1,1,0;AC_NFE=0,0,0;AC_OTH=0,0,0;AC_SAS=0,3,0;AF=3.308e-05,3.308e-05,5.514e-05;AN=90684;AN_AFR=7828;AN_AMR=6546;AN_Adj=79012;AN_EAS=8394;AN_FIN=3354;AN_NFE=39846;AN_OTH=606;AN_SAS=12438;BaseQRankSum=0.736;ClippingRankSum=0.198;DB;DP=1383162;FS=1.848;GQ_HIST=533|4275|91|84|874|13|15|4|2|0|1|0|28889|7105|1334|1168|438|98|75|343,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|2,0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|2,1|2|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0;GQ_MEAN=57.16;GQ_STDDEV=20.29;Het_AFR=0,0,0,0,0,0;Het_AMR=1,0,0,0,0,0;Het_EAS=0,0,0,0,0,0;Het_FIN=0,0,0,0,0,0;Het_NFE=0,0,0,0,0,0;Het_OTH=0,0,0,0,0,0;Het_SAS=0,1,0,0,0,0;Hom_AFR=0,0,0;Hom_AMR=1,0,0;Hom_EAS=0,0,0;Hom_FIN=0,0,0;Hom_NFE=0,0,0;Hom_OTH=0,0,0;Hom_SAS=0,1,0;InbreedingCoeff=0.0345;MQ=30.21;MQ0=0;MQRankSum=-1.231e+00;NCC=25596;QD=9.65;ReadPosRankSum=0.920;VQSLOD=-2.686e+00;culprit=MQ
1	69590	rs141776804	T	A	222918.56	PASS	00_NCBI135_GENEINFO=.;00_NCBI135_VP=050200000004000000000100;AC=110;AC_AFR=1;AC_AMR=102;AC_Adj=103;AC_EAS=0;AC_FIN=0;AC_Het=37;AC_Hom=33;AC_NFE=0;AC_OTH=0;AC_SAS=0;AF=1.172e-03;AN=93836;AN_AFR=7862;AN_AMR=8636;AN_Adj=83862;AN_EAS=8438;AN_FIN=3596;AN_NFE=41696;AN_OTH=630;AN_SAS=13004;BaseQRankSum=0.266;ClippingRankSum=-8.820e-01;DB;DP=1437071;FS=0.000;GQ_HIST=150|4009|30|54|880|37|22|5|6|8|9|15|30709|7422|1387|1208|454|97|77|339,1|2|0|2|0|0|0|0|0|1|0|1|1|0|0|0|0|0|0|68;GQ_MEAN=58.94;GQ_STDDEV=53.31;Het_AFR=1;Het_AMR=36;Het_EAS=0;Het_FIN=0;Het_NFE=0;Het_OTH=0;Het_SAS=0;Hom_AFR=0;Hom_AMR=33;Hom_EAS=0;Hom_FIN=0;Hom_NFE=0;Hom_OTH=0;Hom_SAS=0;InbreedingCoeff=0.0538;MQ=39.62;MQ0=0;MQRankSum=0.559;NCC=23145;QD=18.54;ReadPosRankSum=0.033;VQSLOD=-2.401e-01;culprit=MQ

(...)
```

## History

2018-10-31: add buffered list to speed up things
2017-06-08: more intelligent for AlleleCount.A and AlleleCount.R
2018-07-13: ignore spanning deletions, (for @SolenaSLS)

