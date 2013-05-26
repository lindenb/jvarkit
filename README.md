JVARKIT
=======

Java utilities for Next Generation Sequencing

Pierre Lindenbaum PhD

http://plindenbaum.blogspot.com

@yokofakun	

Dependencies
------------

Tested with java 1.7 and the Picard library 1.91 ( http://sourceforge.net/projects/picard/)


Download & Install
------------------

```bash
git clone "https://github.com/lindenb/jvarkit.git"
cd jvarkit.git
```

edit **build.properties** to configure the project. Something like:

```
picard.jar=/path/to/picard-tools-1.91/picard-1.91.jar
sam.jar=/path/to/picard-tools-1.91/sam-1.91.jar
variant.jar=/path/to/picard-tools-1.91/variant-1.91.jar
tribble.jar=/path/to/picard-tools-1.91/tribble-1.91.jar
```

= Tools =

== Filtering VCF with javascript (rhino) ==
=== Usage ===
```
 java -jar dist/vcffilterjs.jar [option] (vcf|stdin)
```
=== Options ===

* -e (script)

OR

* -f (srcipt-file)


the script binds the following variables:

* **variant** : the current variation;  a org.broadinstitute.variant.variantcontext.VariantContext ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/variantcontext/VariantContext.java )
* **header** : the VCF header org.broadinstitute.variant.vcf.VCFHeader ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/vcf/VCFHeader.java).

the script should return '1' or true if the current VCF file should be printed.

==Example==

the file filter.js

```javascript
/** prints a VARIATION if two samples at least
have a DP<200 */ 
function myfilterFunction()
	{
	var samples=header.genotypeSamples;
	var countOkDp=0;


	for(var i=0; i< samples.size();++i)
		{
		var sampleName=samples.get(i);
		if(! variant.hasGenotype(sampleName)) continue;
		var genotype = variant.genotypes.get(sampleName);
		if( ! genotype.hasDP()) continue;
		var dp= genotype.getDP();
		if(dp < 200 ) countOkDp++;
		}
	return (countOkDp>2)
	}
myfilterFunction();
```

```bash
$ curl -s "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
   java -jar  dist/vcffilterjs.jar  -f filter.js
   
##fileformat=VCFv4.1
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BLANK	NA12878	NA12891	NA12892	NA19238	NA19239	NA19240
chr22	42526449	.	T	A	151.47	.	AC=1;AF=0.071;AN=14;BaseQRankSum=2.662;DP=1226;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=41.2083;MQ=240.47;MQ0=0;MQRankSum=0.578;QD=4.89;ReadPosRankSum=3.611	GT:AD:DP:GQ:PL	0/1:23,8:31:99:190,0,694	0/0:188,0:190:99:0,478,5376	0/0:187,0:187:99:0,493,5322	0/0:247,0:249:99:0,634,6728	0/0:185,0:185:99:0,487,5515	0/0:202,0:202:99:0,520,5857	0/0:181,1:182:99:0,440,5362
chr22	42526634	.	T	C	32.60	.	AC=1;AF=0.071;AN=14;BaseQRankSum=1.147;DP=1225;DS;Dels=0.00;FS=0.000;HRun=0;HaplotypeScore=50.0151;MQ=240.65;MQ0=0;MQRankSum=1.151;QD=1.30;ReadPosRankSum=1.276	GT:AD:DP:GQ:PL	0/1:21,4:25:71:71,0,702	0/0:187,2:189:99:0,481,6080	0/0:233,0:233:99:0,667,7351	0/0:230,0:230:99:0,667,7394	0/0:174,1:175:99:0,446,5469	0/0:194,2:196:99:0,498,6239	0/0:174,0:175:99:0,511,5894
chr22	42527793	rs1080989	C	T	3454.66	.	AC=2;AF=0.167;AN=12;BaseQRankSum=-3.007;DB;DP=1074;DS;Dels=0.01;FS=0.000;HRun=1;HaplotypeScore=75.7865;MQ=209.00;MQ0=0;MQRankSum=3.014;QD=9.36;ReadPosRankSum=0.618	GT:AD:DP:GQ:PL	./.	0/1:72,90:162:99:1699,0,1767	0/1:103,96:202:99:1756,0,2532	0/0:188,0:188:99:0,526,5889	0/0:160,0:160:99:0,457,4983	0/0:197,0:198:99:0,544,6100	0/0:156,0:156:99:0,439,5041

```

=== compile ===

```bash
$ ant vcffilterjs 
Buildfile: /home/lindenb/src/jvarkit-git/build.xml

vcffilterjs:
    (...)

BUILD SUCCESSFUL
Total time: 1 second

```



