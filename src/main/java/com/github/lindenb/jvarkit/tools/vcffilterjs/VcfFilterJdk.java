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
import java.util.List;
import java.util.function.Function;
import com.github.lindenb.jvarkit.lang.InMemoryCompiler;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.samtools.util.IOUtil;
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


The user code is a piece of java code that will be inserted as the body of the following instance of `java.util.function.Function<VariantContext,Object>`:

```java
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.*;
public class VcfFilterJdkCustom1499244002429 implements java.util.function.Function<VariantContext, Object> {
  protected final VCFHeader header;
  public VcfFilterJdkCustom1499244002429(final VCFHeader header) {
  this.header = header;
  }
  @Override
  public Object apply(final VariantContext variant) {
   // user code starts here
   
   // USER CODE is inserted here <-------
   
   // user code ends here
   }
}
```

The program is then compiled in **memory**.

The Function returns an object that can either:


* a boolean : true accept the variant, false reject the variant
* a (VariantContext)[https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html] to replace the current variant
* a (java.util.List<VariantContext>)[https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/variant/variantcontext/VariantContext.html] to replace the current variant with a list of variants.

## See also

* VcfFilterJS

## About Galaxy

At first, this tool is not safe for a public Galaxy server, because the javascript code can access the filesystem.
But you can use the JVM parameter

```
-J-Djava.security.manager
```

to prevent it to access the filesystem. See [http://stackoverflow.com/questions/40177810](http://stackoverflow.com/questions/40177810)


##  Examples


###  Example

filter homozygotes for sample NA12878


```
$ curl -sL "https://raw.github.com/jamescasbon/PyVCF/master/vcf/test/gatk.vcf" |\
	java -jar dist/vcffilterjdk.jar  -e 'return variant.getGenotype("NA12878").isHom();' | grep CHROM -A2

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	BLANK	NA12878	NA12891	NA12892	NA19238	NA19239	NA19240
chr22	42522755	.	C	G	36.98	.	AC=1;AF=0.071;AN=14;BaseQRankSum=-14.866;DP=1527;DS;Dels=0.01;FS=0.000;HRun=0;HaplotypeScore=253.4254;MQ=197.36;MQ0=2;MQRankSum=-10.810;QD=0.15;ReadPosRankSum=-17.244	GT:AD:DP:GQ:PL	0/0:26,1:27:51:0,51,570	0/0:208,40:248:99:0,236,4169	0/0:192,56:249:99:0,114,4292	0/1:179,66:245:75:75,0,3683	0/0:214,32:246:99:0,172,4235	0/0:200,49:249:61:0,61,4049	0/0:195,50:246:32:0,32,3757
chr22	42523003	rs116917064	A	G	7113.55	.	AC=8;AF=0.571;AN=14;BaseQRankSum=6.026;DB;DP=1433;DS;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=101.7894;MQ=182.04;MQ0=0;MQRankSum=-2.501;QD=4.96;ReadPosRankSum=8.294	GT:AD:DP:GQ:PL	0/1:10,2:12:1:1,0,257	1/1:9,173:183:99:2385,273,0	0/1:153,95:249:99:355,0,2355	0/1:140,110:250:99:1334,0,2242	0/1:164,85:249:99:1070,0,2279	0/1:160,90:250:99:1245,0,2300	0/1:156,81:238:99:724,0,2764
```

END_DOC
 */
@Program(
		name="vcffilterjdk",
		description="Filtering VCF with in-memory-compiled java expressions",
		keywords={"vcf","filter","java","jdk"}
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
			
			final  SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header);

			
			final InMemoryCompiler inMemoryCompiler = new InMemoryCompiler();
			
			
			final String javaClassName =VcfFilterJdk.class.getSimpleName()+"Custom"+System.currentTimeMillis();
			
			final StringWriter codeWriter=new StringWriter();
			final PrintWriter pw = new PrintWriter(codeWriter);
			pw.println("import htsjdk.samtools.util.*;");
			pw.println("import htsjdk.variant.variantcontext.*;");
			pw.println("import htsjdk.variant.vcf.*;");
			pw.println("public class "+javaClassName+" implements java.util.function.Function<VariantContext, Object> {");
			pw.println("  protected final VCFHeader header;");
			pw.println("  public "+javaClassName+"(final VCFHeader header) {");
			pw.println("  this.header = header;");
			pw.println("  }");
			pw.println("  @Override");
			pw.println("  public Object apply(final VariantContext variant) {");
			pw.println("   /** user code starts here */");
			pw.println(code);
			pw.println(    "/** user code ends here */");
			pw.println("   }");
			pw.println("}");
			pw.flush();
			
			LOG.debug("#### Compiling :\n" + codeWriter+"\n####");
			
			final Class<?> compiledClass = inMemoryCompiler.compileClass(
					javaClassName,
					codeWriter.toString()
					);
			final Constructor<?> ctor=compiledClass.getDeclaredConstructor(VCFHeader.class);
			@SuppressWarnings("unchecked")
			final Function<VariantContext,Object> filter = (Function<VariantContext,Object>)ctor.newInstance(header);
			
			w.writeHeader(h2);
			while (r.hasNext() && !w.checkError())
				{
				final  VariantContext variation = progress.watch(r.next());

				final Object result = filter.apply(variation);
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
				}
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
