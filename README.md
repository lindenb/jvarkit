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
picard.version=1.91
picard.dir=/home/lindenb/package/picard-tools-${picard.version}
picard.jar=${picard.dir}/picard-${picard.version}.jar
sam.jar=${picard.dir}/sam-${picard.version}.jar
variant.jar=${picard.dir}/variant-${picard.version}.jar
tribble.jar=${picard.dir}/tribble-${picard.version}.jar
```

<h2>Tools</h2>

<h3> Filtering VCF with javascript (rhino) </h3>
<h4>Usage</h4>
```
 java -jar dist/vcffilterjs.jar [option] (vcf|stdin)
```
<h4>Options</h4>

* -e (script)

OR

* -f (srcipt-file)


the script binds the following variables:

* **variant** : the current variation;  a org.broadinstitute.variant.variantcontext.VariantContext ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/variantcontext/VariantContext.java )
* **header** : the VCF header org.broadinstitute.variant.vcf.VCFHeader ( http://sourceforge.net/p/picard/code/HEAD/tree/trunk/src/java/org/broadinstitute/variant/vcf/VCFHeader.java).

the script should return '1' or true if the current VCF file should be printed.

<h4>Example</h4>

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

<h4>compile</h4>

```bash
$ ant vcffilterjs 
Buildfile: /home/lindenb/src/jvarkit-git/build.xml

vcffilterjs:
    (...)

BUILD SUCCESSFUL
Total time: 1 second

```

<h3>SAM4WebLogo</h3>
<h4>Motivation</h4>
"Sequence logo ( http://weblogo.berkeley.edu/logo.cgi ) for different alleles or generated from SAM/BAM" http://www.biostars.org/p/73021

![ScreenShot](https://raw.github.com/lindenb/jvarkit/master/doc/sam2weblogo.png)


<h4>Compilation</h4>

```bash
ant sam4weblogo
```

<h4>Options</h4>

* INPUT=File                        A BAM file to process.  Required. 
* REGION=String Region to observe: chrom:start-end  Required.

<h4>Example</h4>
```bash
$ java -jar dist/sam4weblogo.jar I=path/to/samtools-0.1.18/examples/sorted.bam R=seq1:80-110 2> /dev/null | head -n 50
>B7_593:4:106:316:452/1
TGTTG--------------------------
>B7_593:4:106:316:452a/1
TGTTG--------------------------
>B7_593:4:106:316:452b/1
TGTTG--------------------------
>B7_589:8:113:968:19/2
TGGGG--------------------------
>B7_589:8:113:968:19a/2
TGGGG--------------------------
>B7_589:8:113:968:19b/2
TGGGG--------------------------
>EAS54_65:3:321:311:983/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983a/1
TGTGGG-------------------------
>EAS54_65:3:321:311:983b/1
TGTGGG-------------------------
>B7_591:6:155:12:674/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674a/2
TGTGGGGG-----------------------
>B7_591:6:155:12:674b/2
TGTGGGGG-----------------------
>EAS219_FC30151:7:51:1429:1043/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043a/2
TGTGGGGGGCGCCG-----------------
>EAS219_FC30151:7:51:1429:1043b/2
TGTGGGGGGCGCCG-----------------
>B7_591:5:42:540:501/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501a/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410a/1
TGGGGGGGGCGCAGT----------------
>B7_591:5:42:540:501b/1
TGTGGGGGCCGCAGTG---------------
>EAS192_3:5:223:142:410b/1
TGGGGGGGGCGCAGT----------------
```

<h3>SAM2Tsv</h3>
display a tabular view of each base of the reads vs the reference.
<h4>Compilation</h4>
```bash
ant sam2tsv
```
<h4>Options</h4>
<ul>
<li>IN=File BAM files to process.  This option may be specified 0 or more times.</li>
<li>REGION=String restrict to that region (chr:start-end)  Default value: null. </li>
<li>REF=File Indexed reference  Required. </li>
<li>A={true,false} Use Alignment format.</li>
</ul>
<h4>Example</h4>
```bash
java -jar dist/sam2tsv.jar \
	I=sorted.bam \
	R=genome.fa \
	L=chr1:32944435-32944435


M00491:12:000000000-A3FL3:1:1101:16929:4287	147	1	A	20	chr22	544289	A	M	=
M00491:12:000000000-A3FL3:1:1101:16929:4287	147	2	G	28	chr22	544290	G	M	=
M00491:12:000000000-A3FL3:1:1101:16929:4287	147	3	A	32	chr22	544291	C	M	X
M00491:12:000000000-A3FL3:1:1101:16929:4287	147	4	T	37	chr22	544292	T	M	=
M00491:12:000000000-A3FL3:1:1101:16929:4287	147	5	C	36	chr22	544293	C	M	=
```

<h3>cmpbam: Comparing two or more BAMS</h3>

<h4>Compilation</h4>
```bash
ant cmpbams
```

<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>BAM files to process.  This option must be specified at least 2 times. </td></tr>
<tr><td>REGION=String</td><td>restrict to that region (chr:start-end)  Default value: null. </td></tr>
<tr><td>USESAMFLAG=Boolean</td><td>use SAM Flag when comparing.  Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} </td></tr>
</table>
<h4>Example</h4>
```bash
java -jar dist/cmpbams.ja \
	I=file1.bam \
	I=file2.bam \
	I=file3.bam \
	L=chr1:32944435-32944435
```


<h3>Bam2Raster</h3>
(under development)
<h4>Motivation</h4>
Save a BAM alignment as a PNG image.
<h4>Options</h4>
<table>
<tr><th>Option</th><th>Description</th></tr>
<tr><td>IN=File</td><td>BAM files to process.  Required. </td></tr>
<tr><td>OUT=File</td><td>Image name.  Required. </td></tr>
<tr><td>REGION=String</td><td>restrict to that region (chr:start-end)  Default value: null. </td></tr>
<tr><td>REF=File</td><td>Indexex reference  Required. </td></tr>
<tr><td>WIDTH=Integer</td><td>image width  Default value: 1000. This option can be set to 'null' to clear the default value. </td></tr>
</table>
<h4>Example</h4>
```bash
java -jar dist/bam2raster.jar \
	IN=sorted.bam L=seq1:200-300 \
	OUT=ouput.png \
	R=ex1.fa
```


