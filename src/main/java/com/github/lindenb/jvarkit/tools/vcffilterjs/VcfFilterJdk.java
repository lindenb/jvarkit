/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.vcffilterjs;


import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Constructor;
import java.util.Arrays;
import java.util.Collection;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.OpenJdkCompiler;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.vcf.VariantAttributesRecalculator;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
/**
BEGIN_DOC



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


END_DOC
 */
@Program(
		name="vcffilterjdk",
		description="Filtering VCF with dynamically-compiled java expressions",
		keywords={"vcf","filter","java","jdk"},
		biostars={266201,269854,277820,250212,284083,292710,293314,295902,296145,302217,
				304979,310155,317388,319148,327035,337645,343569,
				347173,351205,351404,354126,302217
				},
		references="\"bioalcidae, samjs and vcffilterjs: object-oriented formatters and filters for bioinformatics files\" . Bioinformatics, 2017. Pierre Lindenbaum & Richard Redon  [https://doi.org/10.1093/bioinformatics/btx734](https://doi.org/10.1093/bioinformatics/btx734).",
		modificationDate="20190222"
		)
public class VcfFilterJdk
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfFilterJdk.class).make();
	@SuppressWarnings("unused")
	private static final Counter<?> _fool_javac=null;
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	
	@Parameter(names={"-F","--filter"},description="If not empty, variants won't be discarded and this name will be used in the FILTER column")
	private String filteredTag = "";
	
	@Parameter(names={"-e","--expression"},description=" (js expression). Optional.")
	private String scriptExpr=null;
	
	@Parameter(names={"-f","--script"},description=" (js file). Optional.")
	private File scriptFile=null;
	
	@Parameter(names={"--nocode"},description=" Don't show the generated code")
	private boolean hideGeneratedCode=false;
	
	@Parameter(names={"--body"},description="user's code is the whole body of the filter class, not just the 'apply' method.")
	private boolean user_code_is_body=false;
	
	@Parameter(names={"--saveCodeInDir"},description="Save the generated java code in the following directory")
	private File saveCodeInDir=null;
	
	@Parameter(names={"-vn","--variable"},description="[20180716] how to name the VariantContext in the code. htsjdk/gatk often use 'vc'.")
	private String variantVariable="variant";
	
	@Parameter(names={"-rc","--recalc"},description="[20180716] recalc attributes like INFO/AF, INFO/AC, INFO/AN... if the number of genotypes has been altered. Recal is not applied if there is no genotype.")
	private boolean recalcAttributes = false;
	
	@Parameter(names={"-xf","--extra-filters"},description="[20180716] extra FILTERs names that will be added in the VCF header and that you can add in the variant using https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContextBuilder.html#filter-java.lang.String- . Multiple separated by space/comma")
	private String extraFilters = "";

			
			
				
				
				
				
				
				public void add(final VariantContext variation) {
					}
				
				
			
			
			
			
	
	public static class AbstractFilter
		extends VcfTools
		implements Function<VariantContext,Object>
		{
		protected final Map<String,Object> userData = new HashMap<>();
		protected final VCFHeader header;
		
		/** Utility: Predicate for filtering the genotypes having at least one ALT allele */
		public final Predicate<Genotype> genotypeHasAltAllele = G->{
			if(G==null) return false;
			final List<Allele> alleles = G.getAlleles();
			for(int i=0;i< alleles.size();i++)
				{
				final Allele a=alleles.get(i);
				if(a.isNoCall()) continue;
				if(a.isReference()) continue;
				return true;
				}
			return false;
			};
		
		protected AbstractFilter(final VCFHeader header) {
			super(header);
			this.header = header;
			}
		@Override
		public Object apply(final VariantContext variant) {
			throw new IllegalStateException("apply(variant) for AbstractFilter is not implemented");
			}
		}
	
	
	public VcfFilterJdk()
		{
		
		}
	
	private String getVariantVariableName() {
		return StringUtil.isBlank(this.variantVariable)?"variant":this.variantVariable;
		}
	
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VCFIterator iter,
			final VariantContextWriter out
			)
		{	
		ProgressFactory.Watcher<VariantContext> progress = null;
		String code = null;
		
		try {
			
			if(this.scriptFile!=null)
				{
				code = IOUtil.slurp(this.scriptFile);
				}
			else
				{
				code = this.scriptExpr;
				}
			final Random rand= new  Random(System.currentTimeMillis());
			final String javaClassName =VcfFilterJdk.class.getSimpleName()+
					"Custom"+ Math.abs(rand.nextInt());
			
			final String generatedClassName= OpenJdkCompiler.getGeneratedAnnotationClassName();
			final StringWriter codeWriter=new StringWriter();
			final PrintWriter pw = new PrintWriter(codeWriter);
			pw.println("import java.util.*;");
			pw.println("import java.util.stream.*;");
			pw.println("import java.util.function.*;");
			pw.println("import htsjdk.samtools.util.*;");
			pw.println("import htsjdk.variant.variantcontext.*;");
			pw.println("import htsjdk.variant.vcf.*;");
	
			if(!StringUtil.isBlank(generatedClassName)) {
				pw.println("@"+generatedClassName+"(value=\""+VcfFilterJdk.class.getSimpleName()+"\",date=\""+ new Iso8601Date(new Date()) +"\")");
				}
			pw.println("public class "+javaClassName+" extends "+AbstractFilter.class.getName().replace('$', '.')+" {");
			pw.println("  public "+javaClassName+"(final VCFHeader header) {");
			pw.println("  super(header);");
			pw.println("  }");
			if(this.user_code_is_body)
				{
				pw.println("   /** user's code starts here */");
				pw.println(code);
				pw.println(    "/** user's code ends here */");
				}
			else
				{
				pw.println("  @Override");
				pw.println("  public Object apply(final VariantContext "+getVariantVariableName()+") {");
				pw.println("   /** user's code starts here */");
				pw.println(code);
				pw.println(    "/** user's code ends here */");
				pw.println("   }");
				}
			pw.println("}");
			pw.flush();
			
			
			if(!this.hideGeneratedCode)
				{
				LOG.debug(" Compiling :\n" + OpenJdkCompiler.beautifyCode(codeWriter.toString()));
				}
			
			if(this.saveCodeInDir!=null)
				{
				PrintWriter cw=null;
				try 
					{
					IOUtil.assertDirectoryIsWritable(this.saveCodeInDir);
					cw = new PrintWriter(new File(this.saveCodeInDir,javaClassName+".java"));
					cw.write(codeWriter.toString());
					cw.flush();
					cw.close();
					cw=null;
					LOG.info("saved "+javaClassName+".java in "+this.saveCodeInDir);
					}
				catch(final Exception err)
					{
					throw new RuntimeIOException(err);
					}
				finally
					{
					CloserUtil.close(cw);
					}
				}
			
			final OpenJdkCompiler compiler = OpenJdkCompiler.getInstance();
			final Class<?> compiledClass = compiler.compileClass(
					javaClassName,
					codeWriter.toString()
					);
			final Constructor<?> constructor = compiledClass.getDeclaredConstructor(VCFHeader.class);
				
				
				
			final VCFHeader header = iter.getHeader();

				
			final AbstractFilter filter_instance;
			;

				
			/* change header */	
			
			
			final VCFHeader h2 = new VCFHeader(header);
			
			/** recalculate INFO/AF, INFO/AN ... variables and 'add' */ 
			final Consumer<VariantContext> recalcAndAdd;
			
			if(this.recalcAttributes && header.hasGenotypingData())
				{
				final VariantAttributesRecalculator recalculator = new VariantAttributesRecalculator();
				recalculator.setHeader(h2);
				recalcAndAdd = V->out.add(recalculator.apply(V));
				}
			else
				{
				recalcAndAdd = V->out.add(V);
				}
			
			final VCFFilterHeaderLine filterHeaderLine = 
				StringUtils.isBlank(filteredTag)?
				null:
				new VCFFilterHeaderLine(this.filteredTag.trim(),"Filtered with "+VcfFilterJdk.class.getSimpleName())
				;
			
			if(filterHeaderLine!=null)
				{
				h2.addMetaDataLine(filterHeaderLine);
				}
			// add genotype filter key if missing
			if(h2.getFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY)==null)
				{
				h2.addMetaDataLine(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY, true));
				}
			
			for(final String xf: Arrays.stream(this.extraFilters.split("[ ,;]+")).
					filter(S->!StringUtil.isBlank(S)).
					collect(Collectors.toSet())
					)
				{
				h2.addMetaDataLine(new VCFFilterHeaderLine(xf,"Custom FILTER inserted with "+VcfFilterJdk.class.getSimpleName()));
				}
			
			try {
				filter_instance = (AbstractFilter)constructor.newInstance(header);
				}
			catch(final Throwable err) {
				LOG.error(err);
				return -1;
				}
			
			/* end change header */
			
			JVarkitVersion.getInstance().addMetaData(this, h2);
			out.writeHeader(h2);
			
			filter_instance.userData.put("first.variant", Boolean.TRUE);
			filter_instance.userData.put("last.variant", Boolean.FALSE);
	
			progress = ProgressFactory.newInstance().dictionary(header).logger(LOG).build();
			while (iter.hasNext() && !out.checkError())
				{				
				final VariantContext variation=progress.apply(iter.next());
				
				/* handle variant */
				final Object result = filter_instance.apply(variation);
				// result is an array of a collection of variants
				if(result!=null && (result.getClass().isArray() || (result instanceof Collection)))
					{
					final  Collection<?> col;
					if(result.getClass().isArray())
						{
						final Object array[]=(Object[])result;
						col= Arrays.asList(array);
						}
					else
						{
						col =( Collection<?>)result;
						}
					// write all of variants
					for(final Object item:col)
						{
						if(item==null) throw new JvarkitException.UserError("item in array is null");
						if(!(item instanceof VariantContext)) throw new JvarkitException.UserError("item in array is not a VariantContext "+item.getClass());
						recalcAndAdd.accept(VariantContext.class.cast(item));
						}
					}
				// result is a VariantContext
				else if(result!=null && (result instanceof VariantContext)) {
					recalcAndAdd.accept(VariantContext.class.cast(result));
					}
				else
					{
					boolean accept=true;
					if(result==null)
						{
						accept=false;
						}
					else if(result instanceof Boolean)
						{
						if(Boolean.FALSE.equals(result)) accept = false;
						}
					else if(result instanceof Number)
						{
						if(((Number)result).intValue()!=1) accept = false;
						}
					else
						{
						LOG.warn("Script returned something that is not a boolean or a number:"+result.getClass());
						accept = false;
						}
					if (!accept)
						{
						if(filterHeaderLine!=null)
							{
							final VariantContextBuilder vcb = new VariantContextBuilder(variation);
							vcb.filter(filterHeaderLine.getID());
							recalcAndAdd.accept(vcb.make());
							}
						continue;
						}
					
					// set PASS filter if needed
					if(filterHeaderLine!=null && !variation.isFiltered())
						{
						recalcAndAdd.accept( new VariantContextBuilder(variation).passFilters().make());
						continue;
						}
					recalcAndAdd.accept(variation);
					}

				
				
				/* end handle variant */
				
				filter_instance.userData.put("first.variant", Boolean.FALSE);
				filter_instance.userData.put("last.variant", !iter.hasNext());
				
				
				final Object stop = filter_instance.userData.get("STOP");
				if(Boolean.TRUE.equals(stop)) break;
				}
			progress.close();
			progress = null;
			out.close();
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			code = null;
			CloserUtil.close(progress);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(this.scriptFile!=null && !StringUtil.isBlank(this.scriptExpr))
			{
			LOG.error("script file and expression both defined");
			return -1;
			}
	
		if(this.scriptFile==null && StringUtil.isBlank(this.scriptExpr))
			{
			LOG.error("script file or expression missing");
			return -1;
			}	
		return doVcfToVcf(args, this.outputFile);
		}

	public static void main(final String[] args) throws Exception
		{
		new VcfFilterJdk().instanceMainWithExit(args);
		}

	}
