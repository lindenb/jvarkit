/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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


import java.io.BufferedWriter;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Constructor;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.OpenJdkCompiler;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.vcf.VcfTools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFRecordCodec;
/**
BEGIN_DOC



## About the script


The user code is a piece of java code that will be inserted as the method apply, or as the body of a com.github.lindenb.jvarkit.tools.vcffilterjs.VcfCustomSortJdk.AbstractComparator class with implements `java.util.Comparator &lt;VariantContext&gt;`:

At the time of writing the documentation, the parent class AbstractFilter is defined as:

```java
public static class AbstractComparator
	extends com.github.lindenb.jvarkit.util.vcf.VcfTools
	implements Comparator<VariantContext>>
	{
	protected final VCFHeader header;
	protected AbstractFilter(final VCFHeader header) {
		super(header);
		this.header = header;
		}
	@Override
	public int compare(final VariantContext vc1 , final VariantContext vc2) {
		throw new IllegalStateException("apply(variant) for AbstractFilter is not implemented");
		}
	// if option --pedigree is defined. Returns an instance of com.github.lindenb.jvarkit.pedigree.Pedigree
	public Pedigree getPedigree();
    public boolean hasPedigree();
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
 9  public class VcfFilterJdkCustom123 extends com.github.lindenb.jvarkit.tools.vcffilterjs.VcfCustomSortJdk.AbstractComparator {
10    public VcfFilterJdkCustom123(final VCFHeader header) {
11    super(header);
12    }
13    @Override
	public int compare(final VariantContext vc1 , final VariantContext vc2) {
15     // user code starts here 
16     user's code is inserted here <===================
17     // user code ends here 
18     }
19  }
```

When the option `--body` is set : the user's code is the whole body (but the constructor) of the class



END_DOC
 */
@Program(
		name="vcfsortjdk",
		description="Sort a  VCF with dynamically-compiled java expressions",
		keywords={"vcf","filter","java","jdk","sort"},
		creationDate="20260426",
		modificationDate="20260429",
		jvarkit_amalgamion =  true,
		menu="VCF Manipulation"
		)
public class VcfCustomSortJdk extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfCustomSortJdk.class);
	@SuppressWarnings("unused")
	private static final Counter<?> _fool_javac=null;
	
	
	@Parameter(names={"--limit"},description="limit output to N variant. -1 == unlimited")
	private long limit_count_variants = -1L;
	
	@Parameter(names={"-e","--expression"},description="The java code expression.")
	private String scriptExpr=null;
	
	@Parameter(names={"-f","--script"},description="The java source code file.")
	private Path scriptPath=null;
	
	@Parameter(names={"--nocode"},description=" Don't show the generated code")
	private boolean hideGeneratedCode=false;
	
	@Parameter(names={"--body"},description="user's code is the whole body of the filter class, not just the 'apply' method.")
	private boolean user_code_is_body=false;
	
	@Parameter(names={"--saveCodeInDir"},description="Save the generated java code in the following directory")
	private Path saveCodeInDir=null;
			
	@Parameter(names={"--tag"},description="Insert the sorting index with this INFO tag. Useful if a standard sort on coordinate is applied after.")
	private String sort_index_tag="IDX";

	
    @ParametersDelegate
    private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();

			
	
	public static class AbstractComparator
		extends VcfTools
		implements Comparator<VariantContext>
		{
		protected final VCFHeader header;
		protected Comparator<VariantContext> defaultComparator;
		
		protected AbstractComparator(final VCFHeader header) {
			super(header);
			this.header = header;
			this.defaultComparator = header.getVCFRecordComparator();
			}
		@Override
		public int compare(final VariantContext vc1,final VariantContext vc2) {
			return this.defaultComparator.compare(vc1, vc2);
			}
		}
	
	
	public VcfCustomSortJdk()
		{
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	@Override
	protected int beforeVcf() {
		if(StringUtils.isBlank(this.sort_index_tag)) {
			LOG.error("INFO sort index tag is blank");
			return -1;
			}
		
		if(this.scriptPath!=null && !StringUtil.isBlank(this.scriptExpr))
			{
			LOG.error("script file and expression both defined");
			return -1;
			}
	
		if(this.scriptPath==null && StringUtil.isBlank(this.scriptExpr))
			{
			LOG.error("script file or expression missing");
			return -1;
			}	
		
		return super.beforeVcf();
		}
	
    @SuppressWarnings("unchecked")
	@Override
    protected int doVcfToVcf(final String inputName, final VCFIterator iter, final VariantContextWriter w) {
    	SortingCollection<VariantContext> sorted=null;
		String code = null;
			try {
				if(this.scriptPath!=null)
					{
					code = IOUtils.slurpPath(this.scriptPath);
					}
				else
					{
					code = this.scriptExpr;
					}
				final Random rand= new  Random(System.currentTimeMillis());
				final String javaClassName =VcfCustomSortJdk.class.getSimpleName()+
						"Custom"+ Math.abs(rand.nextInt());
				
				final StringWriter codeWriter=new StringWriter();
				final PrintWriter pw = new PrintWriter(codeWriter);
				pw.println("import java.util.*;");
				pw.println("import java.util.stream.*;");
				pw.println("import java.util.function.*;");
				pw.println("import htsjdk.samtools.util.*;");
				pw.println("import htsjdk.variant.variantcontext.*;");
				pw.println("import htsjdk.variant.vcf.*;");
		
				pw.println("public class "+javaClassName+" extends "+AbstractComparator.class.getName().replace('$', '.')+" {");
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
					pw.println("  public int compare(final VariantContext vc1,final VariantContext vc2) {");
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
					try(BufferedWriter cw = Files.newBufferedWriter(this.saveCodeInDir.resolve(javaClassName+".java"))) {
						cw.write(codeWriter.toString());
						cw.flush();
						}
					catch(final Exception err)
						{
						throw new RuntimeIOException(err);
						}
					LOG.info("saved "+javaClassName+".java in "+this.saveCodeInDir);
					}
				
				final OpenJdkCompiler compiler = OpenJdkCompiler.getInstance();
				final Class<?> compiledClass = compiler.compileClass(
						javaClassName,
						codeWriter.toString()
						);
				final Constructor<?> constructor = compiledClass.getDeclaredConstructor(VCFHeader.class);
					
					
			
			
			final VCFHeader header= iter.getHeader();
		
				
			final Comparator<VariantContext> cmp;
			try {
				cmp = (Comparator<VariantContext>)constructor.newInstance(header);
				}
			catch(final Throwable err) {
				LOG.error(err);
				return -1;
				}
			final VCFHeader header2 = new VCFHeader(header);
			final VCFInfoHeaderLine info_idx= new VCFInfoHeaderLine(this.sort_index_tag, 1, VCFHeaderLineType.Integer,"Order after "+getProgramName());
			header2.addMetaDataLine(info_idx);
			JVarkitVersion.getInstance().addMetaData(getClass().getSimpleName(), header2);
	
		
			sorted=SortingCollection.newInstance(
					VariantContext.class,
					new VCFRecordCodec(header),
					cmp,
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sorted.setDestructiveIteration(true);
			while(iter.hasNext())
				{
				sorted.add(iter.next());
				}
			
			sorted.doneAdding();
			
			w.writeHeader(header2);
			try(CloseableIterator<VariantContext> iter2 =sorted.iterator()) {
				long idx=1L;
				while(iter2.hasNext())
					{
					if(limit_count_variants>=0L && idx> limit_count_variants) break;
					w.add(new VariantContextBuilder(iter2.next()).attribute(info_idx.getID(), idx).make());
					++idx;
					}
				}
			sorted.doneAdding();
			sorted=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			try {
				if(sorted!=null) sorted.cleanup();
				} catch(Exception err){}
			}
		}
	
	

	public static void main(final String[] args) throws Exception
		{
		new VcfCustomSortJdk().instanceMainWithExit(args);
		}

	}
