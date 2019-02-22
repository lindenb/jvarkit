# VcfFilterJdk

![Last commit](https://img.shields.io/github/last-commit/lindenb/jvarkit.png)

Filtering VCF with dynamically-compiled java expressions


## Usage

```
Usage: vcffilterjdk [options] Files
  Options:
    --body
      user's code is the whole body of the filter class, not just the 'apply' 
      method. 
      Default: false
    -e, --expression
       (js expression). Optional.
    -xf, --extra-filters
      [20180716] extra FILTERs names that will be added in the VCF header and 
      that you can add in the variant using https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContextBuilder.html#filter-java.lang.String- 
      . Multiple separated by space/comma
      Default: <empty string>
    -F, --filter
      If not empty, variants won't be discarded and this name will be used in 
      the FILTER column
      Default: <empty string>
    -h, --help
      print help and exit
    --helpFormat
      What kind of help. One of [usage,markdown,xml].
    --nocode
       Don't show the generated code
      Default: false
    -o, --output
      Output file. Optional . Default: stdout
    -rc, --recalc
      [20180716] recalc attributes like INFO/AF, INFO/AC, INFO/AN... if the 
      number of genotypes has been altered. Recal is not applied if there is 
      no genotype.
      Default: false
    --saveCodeInDir
      Save the generated java code in the following directory
    -f, --script
       (js file). Optional.
    -vn, --variable
      [20180716] how to name the VariantContext in the code. htsjdk/gatk often 
      use 'vc'.
      Default: variant
    --version
      print version and exit

```


## Keywords

 * vcf
 * filter
 * java
 * jdk



## See also in Biostars

 * [https://www.biostars.org/p/266201](https://www.biostars.org/p/266201)
 * [https://www.biostars.org/p/269854](https://www.biostars.org/p/269854)
 * [https://www.biostars.org/p/277820](https://www.biostars.org/p/277820)
 * [https://www.biostars.org/p/250212](https://www.biostars.org/p/250212)
 * [https://www.biostars.org/p/284083](https://www.biostars.org/p/284083)
 * [https://www.biostars.org/p/292710](https://www.biostars.org/p/292710)
 * [https://www.biostars.org/p/293314](https://www.biostars.org/p/293314)
 * [https://www.biostars.org/p/295902](https://www.biostars.org/p/295902)
 * [https://www.biostars.org/p/296145](https://www.biostars.org/p/296145)
 * [https://www.biostars.org/p/302217](https://www.biostars.org/p/302217)
 * [https://www.biostars.org/p/304979](https://www.biostars.org/p/304979)
 * [https://www.biostars.org/p/310155](https://www.biostars.org/p/310155)
 * [https://www.biostars.org/p/317388](https://www.biostars.org/p/317388)
 * [https://www.biostars.org/p/319148](https://www.biostars.org/p/319148)
 * [https://www.biostars.org/p/327035](https://www.biostars.org/p/327035)
 * [https://www.biostars.org/p/337645](https://www.biostars.org/p/337645)
 * [https://www.biostars.org/p/343569](https://www.biostars.org/p/343569)
 * [https://www.biostars.org/p/347173](https://www.biostars.org/p/347173)
 * [https://www.biostars.org/p/351205](https://www.biostars.org/p/351205)
 * [https://www.biostars.org/p/351404](https://www.biostars.org/p/351404)
 * [https://www.biostars.org/p/354126](https://www.biostars.org/p/354126)
 * [https://www.biostars.org/p/302217](https://www.biostars.org/p/302217)


## Compilation

### Requirements / Dependencies

* java [compiler SDK 11](https://jdk.java.net/11/). Please check that this java is in the `${PATH}`. Setting JAVA_HOME is not enough : (e.g: https://github.com/lindenb/jvarkit/issues/23 )


### Download and Compile

```bash
$ git clone "https://github.com/lindenb/jvarkit.git"
$ cd jvarkit
$ ./gradlew vcffilterjdk
```

The java jar file will be installed in the `dist` directory.

## Source code 

[https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfFilterJdk.java](https://github.com/lindenb/jvarkit/tree/master/src/main/java/com/github/lindenb/jvarkit/tools/vcffilterjs/VcfFilterJdk.java)


## Contribute

- Issue Tracker: [http://github.com/lindenb/jvarkit/issues](http://github.com/lindenb/jvarkit/issues)
- Source Code: [http://github.com/lindenb/jvarkit](http://github.com/lindenb/jvarkit)

## License

The project is licensed under the MIT license.

## Citing

Should you cite **vcffilterjdk** ? [https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md](https://github.com/mr-c/shouldacite/blob/master/should-I-cite-this-software.md)

The current reference is:

 * "bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734).




## About the script


The user code is a piece of java code that will be inserted as the method apply, or as the body of a com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk.AbstractFilter class with implements `java.util.function.Function<VariantContext,Object>`:

At the time of writing the documentation, the parent class AbstractFilter is defined as:

```java
public static class AbstractFilter
	extends com.github.lindenb.jvarkit.util.vcf.VcfTools
	implements Function<VariantContext,Object>
	{
	protected final Map<String,Object> userData = new HashMap<>();
	protected final VCFHeader header;
	protected AbstractFilter(final VCFHeader header) {
		super(header);
		this.header = header;
		}
	@Override
	public Object apply(final VariantContext variant) {
		throw new IllegalStateException("apply(variant) for AbstractFilter is not implemented");
		}
	}
```

where

* the base class VcfTools contains some utilities for parsing VEP/SNPEFF annotations and detecting Mendelian Violations. see [https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java](https://github.com/lindenb/jvarkit/blob/master/src/main/java/com/github/lindenb/jvarkit/util/vcf/VcfTools.java).
* 'header' is the VCF header
* 'userData' is a placeHolder where the user is free to put things.

'userData' will be filled with the following properties:

* `<"first.variant",Boolean>` the current variant is the first variant in the VCF
* `<"last.variant",Boolean>` the current variant is the last variant in the VCF

If the user puts `<"STOP",Boolean.TRUE>` in `userData` the scanning of the VCF will be aborted without error.

The user code will be inserted in the following java code:


```
 1  import java.util.*;
 2  import java.util.stream.*;
 3  import java.util.function.*;
 4  import htsjdk.samtools.util.*;
 5  import htsjdk.variant.variantcontext.*;
 6  import htsjdk.variant.vcf.*;
 7  import  javax.annotation.processing.Generated;
 8  @Generated("VcfFilterJdk")
 9  public class VcfFilterJdkCustom123 extends com.github.lindenb.jvarkit.tools.vcffilterjs.VcfFilterJdk.AbstractFilter {
10    public VcfFilterJdkCustom123(final VCFHeader header) {
11    super(header);
12    }
13    @Override
14    public Object apply(final VariantContext variant) {
15     // user code starts here 
16     user's code is inserted here <===================
17     // user code ends here 
18     }
19  }
```

When the option `--body` is set : the user's code is the whole body (but the constructor) of the class



The program is then compiled in **memory**.

The method `apply` returns an object that can either:

* a boolean : true accept the variant, false reject the variant
* a [VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) to replace the current variant
* a [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html)<[VariantContext](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html) > to replace the current variant with a list of variants.

## See also

* VcfFilterJS . Slower, using javascript syntax (rhino engine)

## About Galaxy

At first, this tool is not safe for a public Galaxy server, because the javascript code can access the filesystem.
But you can use the JVM parameter

```
-J-Djava.security.manager
```

to prevent it to access the filesystem. See [http://stackoverflow.com/questions/40177810](http://stackoverflow.com/questions/40177810)


##  Examples


###  Example

>  I'm interested in finding all sites (regardless of genotype call heterozygous or homozygous) where at least one of the alternative alleles have an AD value (Allelic Depth) greater than 10,

see [https://bioinformatics.stackexchange.com/questions/974/](https://bioinformatics.stackexchange.com/questions/974/)

```
java -jar dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.hasAD() && java.util.Arrays.stream(G.getAD()).skip(1).filter(AD->AD>10).findAny().isPresent()).findAny().isPresent();' 
```

###  Example

filter homozygotes for sample NA12878


```
$ curl -sL "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
	java -jar dist/vcffilterjdk.jar  -e 'return variant.getGenotype("NA12878").isHom();' | grep CHROM -A2

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BLANK	NA12878	NA12891	NA12892	NA19238	NA19239	NA19240
chr22	42522755	.	C	G	36.98	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-14.866;DP=1527;DS;Dels=0.01;FS=0.000;HRun=0;HaplotypeScore=253.4254;MQ=197.36;MQ0=2;MQRankSum=-10.810;QD=0.15;ReadPosRankSum=-17.244	GT:AD:DP:GQ:PL	0/0:26,1:27:51:0,51,570	0/0:208,40:248:99:0,236,4169	0/0:192,56:249:99:0,114,4292	0/1:179,66:245:75:75,0,3683	0/0:214,32:246:99:0,172,4235	0/0:200,49:249:61:0,61,4049	0/0:195,50:246:32:0,32,3757
chr22	42523003	rs116917064	A	G	7113.55	.	AC=8;AF=0.571;AN=14;BaseQRankSum=6.026;DB;DP=1433;DS;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=101.7894;MQ=182.04;MQ0=0;MQRankSum=-2.501;QD=4.96;ReadPosRankSum=8.294	GT:AD:DP:GQ:PL	0/1:10,2:12:1:1,0,257	1/1:9,173:183:99:2385,273,0	0/1:153,95:249:99:355,0,2355	0/1:140,110:250:99:1334,0,2242	0/1:164,85:249:99:1070,0,2279	0/1:160,90:250:99:1245,0,2300	0/1:156,81:238:99:724,0,2764
```

### Example

first and second genotype are not the same:

```
java -jar dist/vcffilterjdk.jar -e 'return !variant.getGenotype(0).sameGenotype(variant.getGenotype(1));' 
```

### Example

at least 3 samples have a DP greater than 30

```
java -jar dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.hasDP() && G.getDP()>30).limit(3).count()> 2;' 
```

### Example

Variant is annotated with SO:0001818 or its children ( protein_altering_variant )

```
$ java -jar dist/vcffilterjdk.jar -e 'return this.hasSequenceOntologyLabel(variant,"protein_altering_variant");' 
```

or

```
java -jar dist/vcffilterjdk.jar -e 'return this.hasSequenceOntologyAccession(variant,"SO:0001818");' 
```

### Example

Unphase a VCF file

```
java -jar dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes(variant.getGenotypes().stream().map(G->new GenotypeBuilder(G).phased(false).make()).collect(Collectors.toList())).make();' input.vcf
```

## Example

Change haploid to diploid

note: things like 'AF' are not fixed.

```
$ wget -q -O - "https://raw.githubusercontent.com/CostaLab/practical_SS2015/598ea0dddf2ef073a55ae21bc6d39ac2172eb617/data_analysis/organisms/escherichia_coli/O157H7_Sakai/IonTorrentPGM_mem/SRX185723/SRR566635/SRR566635-snps.vcf" | java -jar dist/vcffilterjdk.jar -e 'return new VariantContextBuilder(variant).genotypes(variant.getGenotypes().stream().map(G->!G.isCalled()?GenotypeBuilder.createMissing(G.getSampleName(),2):G).map(G->G.isCalled() && G.getPloidy()==1?new GenotypeBuilder(G).alleles(Arrays.asList(G.getAllele(0),G.getAllele(0))).make():G).collect(Collectors.toList())).attribute("AC",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 && !G.getAllele(0).isReference()?2:G.getAlleles().size()).sum()).attribute("AN",variant.getGenotypes().stream().mapToInt(G->G.isCalled() && G.getPloidy()==1 ?2:G.getAlleles().size()).sum()).make();'
```

## Example

set DP missing in genotypes when only AD is present:

```
curl -s "http://www.icbi.at/svnsimplex/CommonBioCommands/tags/simplex-1.0/CommonBioCommands/testdata/vcf/AMS1_converted_filtered_short_chr1.vcf" |\
awk '/^#CHROM/ {printf("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"\">\n");} {print;}' | \
java -jar dist/vcffilterjdk.jar  -e 'List<Genotype> gl = variant.getGenotypes().stream().map(G->{int ad[]=G.getAD();if(ad==null || ad.length==0) return G; int c= Arrays.stream(ad).sum(); return new GenotypeBuilder(G).DP(c).make();}).collect(Collectors.toList());return new VariantContextBuilder(variant).genotypes(gl).make();'
```

## Example

> How can I access the n-th item in the m-th sample"

```
$ java -jar dist/vcffilterjdk.jar -e 'Genotype G=variant.getGenotype(0); return G.hasAD() && G.getAD().length>1 &&  G.getAD()[1]>3;' input.vcf
```

## Example

> Identifying variants differing between control/treatment

```
$ java -jar dist/vcffilterjdk.jar --body -e 'List<String> cases = null,controls=null;  public Object apply(final VariantContext variant) { if(cases==null) try {cases=IOUtil.slurpLines(new java.io.File("treat.txt")) ; controls= IOUtil.slurpLines(new java.io.File("control.txt")); } catch(Exception e) {throw new RuntimeIOException(e);}; for(final String S1:cases) {final Genotype G1=variant.getGenotype(S1); if(G1==null ||!G1.isCalled()) continue;for(final String S2:controls) {final Genotype G2=variant.getGenotype(S2); if(G2==null ||!G2.isCalled()) continue;   if(G1.sameGenotype(G2)) return false;}} return true;}' input.vcf```
```

## Example

> Retain sites where atleast 80% of the individuals had at least depth DP >= 10 and GQ>=20 irrespective of the reference or non-reference allele

```
java -jar dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().filter(G->G.getDP()>=10 && G.getGQ()>=20).count()/(double)variant.getNSamples() > 0.8;' input.vcf
```

> Retain sites where at least one sample has the non-reference allele with DP>= 10 and GQ >= 20.

```
java -jar dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().anyMatch(G->G.getDP()>=10 && G.getGQ()>=20 && G.getAlleles().stream().anyMatch(A->A.isCalled() && !A.isReference())) ;'
```

## Example

creating a simple GUI for **vcffilterjdk** using [zenity](https://en.wikipedia.org/wiki/Zenity)

see [https://gist.github.com/lindenb/4465c0e822b175f3428029526beef80c](https://gist.github.com/lindenb/4465c0e822b175f3428029526beef80c) , [https://www.biostars.org/p/296145/](https://www.biostars.org/p/296145/)

![capture](https://gist.githubusercontent.com/lindenb/4465c0e822b175f3428029526beef80c/raw/3510261585a1fc8858c2fb54caba2d1c43d72918/Screenshot_A.png)


## Example

updating AF and MAF fields:

```
$ gunzip -c  input.vcf.gz |\
 awk '/^#CHROM/ {printf("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Min Allele Frequency\">\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");} {print}' |\
 java -jar dist/vcffilterjdk.jar -e 'VariantContextBuilder vcb = new VariantContextBuilder(variant); float ac = variant.getAttributeAsInt("AN",0); if(ac>0) { List<Float> af = variant.getAttributeAsIntList("AC",0).stream().map(N->N/ac).collect(Collectors.toList());vcb.attribute("AF",af);vcb.attribute("MAF",af.stream().mapToDouble(X->X.floatValue()).min().orElse(-1.0) );} return vcb.make();'
```

## Example

Set attribute "AA=at|..."  to upper case in vcf file

```
java -jar dist/vcffilterjdk.jar -e 'if(!variant.hasAttribute("AA")) return variant; String AA= variant.getAttributeAsString("AA",""); int pipe=AA.indexOf("|"); AA= AA.substring(0,pipe).toUpperCase()+AA.substring(pipe); return new VariantContextBuilder(variant).attribute("AA",AA).make();'
```

## Example

select LUMPY-SV familial structural variations


run with the option '--body'

```java
private final List<Set<String>> sampleSetList = Arrays.asList(
		new HashSet<>(Arrays.asList("S1","S2")),
		new HashSet<>(Arrays.asList("S3","S4","S5")),
		new HashSet<>(Arrays.asList("S6","S7","S8","S9"))
		);

private boolean isSV(final Genotype g) {
	if(g.getAttributeAsInt("SU",0)>0) return true;
	if(g.getAttributeAsInt("SR",0)>0) return true;
	return false;	
	}


private boolean validateSet(final VariantContext V,final Set<String> sampleSet) {
	if(sampleSet.stream().
		map(N->V.getGenotype(N)).
		anyMatch(G->!isSV(G)) )
		{
		return false;
		}

	if(V.getGenotypes().stream().
		filter(G->!sampleSet.contains(G.getSampleName())).
		anyMatch(G->isSV(G)))
		{
		return false;
		}
	return true;
	}


@Override
public Object apply(final VariantContext V) {
	for(final Set<String> sampleSet : this.sampleSetList)
		{
		if(validateSet(V,sampleSet)) return true;
		}
	return false;
	}
```

## Example

> What I want is to assign ./. missing genotypes for sample-level genotypes in a VCF, if they fail to pass defined AD ratio (so, keeping the variant itself still if it is present in at least one sample with desired AD ratio thresholds). 

```java
return new VariantContextBuilder(variant).
    genotypes(
        variant.
        getGenotypes().
        stream().
        map( G->{
            if(G.hasAD()) {
                final int A[]= G.getAD();
                if(A.length>1 &&  A[0]>0 && A[1]/(double)A[0]> 0.33) return G;
                }
            return  GenotypeBuilder.createMissing(G.getSampleName(),2);
            }).
        collect(Collectors.toList())
        ).
    make();
```

## Example

> I've filtered variants based on quality and depth using bcftools, but I also want to filter them that were called within 5 base pairs of each other. Is there any tools that can do this?

```
$ java -jar dist/vcffilterjdk.jar --body -e 'String prev_contig=null;int prev_end=-1; public Object apply(final VariantContext vc) {boolean ret=!(vc.getContig().equals(prev_contig)  && vc.getStart() - prev_end <=5); prev_contig=vc.getContig();prev_end=vc.getEnd();return ret; }'  in.vcf
```

## Example

> how to filter a multi-VCF for sex-specific genotypes?

see https://bioinformatics.stackexchange.com/questions/5518

```
java -jar dist/vcffilterjdk.jar -e 'return variant.getGenotypes().stream().allMatch(G->(G.getSampleName().endsWith("M") && G.isHet()) || (G.getSampleName().endsWith("F") && G.isHomRef())); ' input.vcf
```

## History:

  * 20190222 : removed some jaxb stuff
  * 201901 : OpenJdk doesn't support anymore in-memory compiling. Switching to OpenJdkCompiler


