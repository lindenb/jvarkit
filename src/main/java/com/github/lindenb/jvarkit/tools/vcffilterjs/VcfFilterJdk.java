/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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
import java.util.function.Function;

import com.github.lindenb.jvarkit.lang.InMemoryCompiler;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
 7  import javax.annotation.Generated;
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


END_DOC
 */
@Program(
		name="vcffilterjdk",
		description="Filtering VCF with in-memory-compiled java expressions",
		keywords={"vcf","filter","java","jdk"},
		biostars={266201}
		)
public class VcfFilterJdk
	extends Launcher
	{	
	private static final Logger LOG = Logger.build(VcfFilterJdk.class).make();
	
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
	
	
	public static class AbstractFilter
		extends VcfTools
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
	
	
	public VcfFilterJdk()
		{
		
		}
	
	@Override
	protected int doVcfToVcf(final String inputName,final VcfIterator r,final VariantContextWriter w) {
		try
			{
			final String code;
			
			if(this.scriptFile!=null)
				{
				code = IOUtil.slurp(this.scriptFile);
				}
			else
				{
				code = this.scriptExpr;
				}
			
			final VCFHeader header = r.getHeader();
			final VCFHeader h2 = new VCFHeader(header);
			addMetaData(h2);
			
			final VCFFilterHeaderLine filterHeaderLine = (filteredTag.trim().isEmpty()?null:
				new VCFFilterHeaderLine(this.filteredTag.trim(),"Filtered with "+getProgramName())
				);
			
			
			if(filterHeaderLine!=null) h2.addMetaDataLine(filterHeaderLine);
			

			
			
			final Random rand= new  Random(System.currentTimeMillis());
			final String javaClassName =VcfFilterJdk.class.getSimpleName()+
					"Custom"+ Math.abs(rand.nextInt());
			
			final StringWriter codeWriter=new StringWriter();
			final PrintWriter pw = new PrintWriter(codeWriter);
			pw.println("import java.util.*;");
			pw.println("import java.util.stream.*;");
			pw.println("import java.util.function.*;");
			pw.println("import htsjdk.samtools.util.*;");
			pw.println("import htsjdk.variant.variantcontext.*;");
			pw.println("import htsjdk.variant.vcf.*;");
			pw.println("import javax.annotation.Generated;");

			pw.println("@Generated(value=\""+VcfFilterJdk.class.getSimpleName()+"\",date=\""+ new Iso8601Date(new Date()) +"\")");
			pw.println("public class "+javaClassName+" extends "+AbstractFilter.class.getName().replace('$', '.')+" {");
			pw.println("  public "+javaClassName+"(final VCFHeader header) {");
			pw.println("  super(header);");
			pw.println("  }");
			if(user_code_is_body)
				{
				pw.println("   /** user's code starts here */");
				pw.println(code);
				pw.println(    "/** user's code ends here */");
				}
			else
				{
				pw.println("  @Override");
				pw.println("  public Object apply(final VariantContext variant) {");
				pw.println("   /** user's code starts here */");
				pw.println(code);
				pw.println(    "/** user's code ends here */");
				pw.println("   }");
				}
			pw.println("}");
			pw.flush();
			
			
			if(!hideGeneratedCode)
				{
				LOG.debug(" Compiling :\n" + InMemoryCompiler.beautifyCode(codeWriter.toString()));
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
					LOG.error(err);
					return -1;
					}
				finally
					{
					CloserUtil.close(cw);
					}
				}
			
			final InMemoryCompiler inMemoryCompiler = new InMemoryCompiler();
			final Class<?> compiledClass = inMemoryCompiler.compileClass(
					javaClassName,
					codeWriter.toString()
					);
			final Constructor<?> ctor=compiledClass.getDeclaredConstructor(VCFHeader.class);
			final AbstractFilter filter = (AbstractFilter)ctor.newInstance(header);
			
			w.writeHeader(h2);
			
			filter.userData.put("first.variant", Boolean.TRUE);
			filter.userData.put("last.variant", Boolean.FALSE);

			final  SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
			while (r.hasNext() && !w.checkError())
				{
				final  VariantContext variation = progress.watch(r.next());

				final Object result = filter.apply(variation);
				
				filter.userData.put("first.variant", Boolean.FALSE);
				filter.userData.put("last.variant", !r.hasNext());
				
				
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
						w.add(VariantContext.class.cast(item));
						}
					}
				// result is a VariantContext
				else if(result!=null && (result instanceof VariantContext)) {
					w.add(VariantContext.class.cast(result));
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
							w.add(vcb.make());
							}
						continue;
						}
					
					// set PASS filter if needed
					if(filterHeaderLine!=null && !variation.isFiltered())
						{
						w.add( new VariantContextBuilder(variation).passFilters().make());
						continue;
						}
					
					w.add(variation);
					}
				
				final Object stop = filter.userData.get("STOP");
				if(Boolean.TRUE.equals(stop)) break;
				
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{

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
		
		try 
			{			
			return doVcfToVcf(args, this.outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(String[] args) throws Exception
		{
		new VcfFilterJdk().instanceMainWithExit(args);
		}

	}
