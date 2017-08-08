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
package com.github.lindenb.jvarkit.tools.samjs;

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

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.InMemoryCompiler;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;


/**
BEGIN_DOC

## Motivation

Filters a BAM using a java expression compiled in memory.


## About the script


The user code is a piece of java code that will be inserted as the method test, or as the body of a com.github.lindenb.jvarkit.tools.samjs.SamJdk.AbstractFilter class with implements `java.util.function.Predicate<SAMRecord>`:

At the time of writing the documentation, the parent class AbstractFilter is defined as:

```java
public static class AbstractFilter
		implements Function<SAMRecord,Object>
		{
		// hashmap, the user is free to use It 
		protected final Map<String,Object> userData = new HashMap<>();
		// input SAM header
		protected final SAMFileHeader header;
		protected AbstractFilter(final SAMFileHeader header) {
			this.header = header;
			}
		@Override
		public Object apply(final SAMRecord record) {
			throw new IllegalStateException("apply(record) for AbstractFilter is not implemented");
			}
		}
```

where

* 'header' is the input SAM header
* 'userData' is a placeHolder where the user is free to put things.

The user code will be inserted in the following java code:


```
 1  import java.util.*;
 2  import java.util.stream.*;
 3  import java.util.function.*;
 4  import htsjdk.samtools.*;
 5  import htsjdk.samtools.util.*;
 6  import javax.annotation.Generated;
 7  @Generated(value="SamJdk",date="2017-08-07T14:48:39+0200")
 8  public class SamJdkCustom756098808 extends com.github.lindenb.jvarkit.tools.samjs.SamJdk.AbstractFilter {
 9    public SamJdkCustom756098808(final SAMFileHeader header) {
10    super(header);
11    }
12    @Override
13    public boolean test(final SAMRecord record) {
14     // user's code starts here 
15     return record.getContig()==null;
16     //user's code ends here 
17     }
18  }
```


When the option `--body` is set : the user's code is the whole body (but the constructor) of the class


The program is then compiled in **memory**.

The method `apply` returns an object that can either:

* a boolean : true accept the read, false reject the record
* a [SAMRecord](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html) to replace the current read
* a [java.util.List](https://docs.oracle.com/javase/8/docs/api/java/util/List.html)<[SAMRecord](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecord.html) > to replace the current record with a list of records.


## Example

### Example 1


get a SAM where the  read OR the mate is unmapped

```bash
java -jar dist/samjdk.jar  \
	-e "return record.getReadUnmappedFlag() || record.getMateUnmappedFlag();" \
	ex1.bam

@HD	VN:1.4	SO:unsorted
@SQ	SN:seq1	LN:1575
@SQ	SN:seq2	LN:1584
B7_591:4:96:693:509	73	seq1	1	99	36M	*	0	0	CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCG	<<<<<<<<<<<<<<<;<<<<<<<<<5<<<<<;:<;7	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
EAS54_65:7:152:368:113	73	seq1	3	99	35M	*	0	0	CTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGT	<<<<<<<<<<0<<<<655<<7<<<:9<<3/:<6):H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
EAS51_64:8:5:734:57	137	seq1	5	99	35M	*	0	0	AGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCC	<<<<<<<<<<<7;71<<;<;;<7;<<3;);3*8/5H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:66
B7_591:1:289:587:906	137	seq1	6	63	36M	*	0	0	GTGGCTCATTGTAATTTTTTGTTTTAACTCTTCTCT	(-&----,----)-)-),'--)---',+-,),''*,	H0:i:0	H1:i:0	MF:i:130	NM:i:5	UQ:i:38	Aq:i:63
EAS56_59:8:38:671:758	137	seq1	9	99	35M	*	0	0	GCTCATTGTAAATGTGTGGTTTAACTCGTCCATGG	<<<<<<<<<<<<<<<;<;7<<<<<<<<7<<;:<5%H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:72
EAS56_61:6:18:467:281	73	seq1	13	99	35M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCCTGGCCCA	<<<<<<<<;<<<8<<<<<;8:;6/686&;(16666H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:5	Aq:i:39
EAS114_28:5:296:340:699	137	seq1	13	99	36M	*	0	0	ATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAG	<<<<<;<<<;<;<<<<<<<<<<<8<8<3<8;<;<0;	H0:i:1	H1:i:0	MF:i:18	NM:i:0	UQ:i:0	Aq:i:73
B7_597:6:194:894:408	73	seq1	15	99	35M	*	0	0	TGTAAATGTGTGGTTTAACTCGTCCATTGCCCAGC	<<<<<<<<<7<<;<<<<;<<<7;;<<<*,;;572<H0:i:0	H1:i:1	MF:i:18	NM:i:1	UQ:i:9	Aq:i:43
EAS188_4:8:12:628:973	89	seq1	18	75	35M	*	0	0	AAATGTGTGGTTTAACTCGTCCATGGCCCAGCATT	==;=:;:;;:====;=;===:=======;==;===H0:i:1	H1:i:0	MF:i:64	NM:i:0	UQ:i:0	Aq:i:0
(...)
```

### Example 2

remove reads with indels:

```
java -jar dist/samjdk.jar -e 'if(record.getReadUnmappedFlag()) return false; Cigar cigar=record.getCigar();if(cigar==null) return false; for(int i=0;i< cigar.numCigarElements();++i) {if(cigar.getCigarElement(i).getOperator().isIndelOrSkippedRegion()) return false; } return true;' input.bam
```


END_DOC
*/
@Program(name="samjdk",
	description="Filters a BAM using a java expression compiled in memory.",
	keywords={"sam","bam","java","jdk","filter"},
	biostars={}
	)
public class SamJdk
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SamJdk.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-X","--fail"},description="Save dicarded reads in that file")
	private File failingReadsFile = null;

	@Parameter(names={"-N","--limit"},description="limit to 'N' records (for debugging).")
	private long LIMIT = -1L ;

	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	@Parameter(names={"-e","--expression"},description="javascript expression")
	private String scriptExpr=null;
	@Parameter(names={"-f","--file"},description="javascript file")
	private File scriptFile =null;
	private SAMFileWriter failingReadsWriter=null;
	
	@Parameter(names={"--nocode"},description=" Don't show the generated code")
	private boolean hideGeneratedCode=false;
	@Parameter(names={"--body"},description="user's code is the whole body of the filter class, not just the 'apply' method.")
	private boolean user_code_is_body=false;
	@Parameter(names={"--saveCodeInDir"},description="Save the generated java code in the following directory")
	private File saveCodeInDir=null;


	public static class AbstractFilter
			implements Function<SAMRecord,Object>
			{
			/** hashmap, the user is free to use It */
			protected final Map<String,Object> userData = new HashMap<>();
			/** input SAM header */
			protected final SAMFileHeader header;
			protected AbstractFilter(final SAMFileHeader header) {
				this.header = header;
				}
			@Override
			public Object apply(final SAMRecord record) {
				throw new IllegalStateException("apply(record) for AbstractFilter is not implemented");
				}
			}

	
	
	
	public SamJdk()
		{
		}

	/* open failing bam if it was not already open */
	private void openFailing(final SAMFileHeader h)
		{
		if(this.failingReadsFile==null) return;
		if(this.failingReadsWriter==null)
			{
			LOG.info("Writing failings to "+ this.failingReadsFile);
			final SAMFileHeader h2= h.clone();
			this.failingReadsWriter=this.writingBamArgs.openSAMFileWriter(failingReadsFile, h2,true);
			}
		}

	private void failing(final SAMRecord rec,final SAMFileHeader h)
		{
		openFailing(h);
		if(failingReadsWriter!=null) failingReadsWriter.addAlignment(rec);
		}

	@Override
	public int doWork(final List<String> args) {
		SAMRecordIterator iter=null;
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
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

			final Random rand= new  Random(System.currentTimeMillis());
			final String javaClassName =SamJdk.class.getSimpleName()+
					"Custom"+ Math.abs(rand.nextInt());
			
			final StringWriter codeWriter=new StringWriter();
			final PrintWriter pw = new PrintWriter(codeWriter);
			pw.println("import java.util.*;");
			pw.println("import java.util.stream.*;");
			pw.println("import java.util.function.*;");
			pw.println("import htsjdk.samtools.*;");
			pw.println("import htsjdk.samtools.util.*;");
			pw.println("import javax.annotation.Generated;");

			pw.println("@Generated(value=\""+SamJdk.class.getSimpleName()+"\",date=\""+ new Iso8601Date(new Date()) +"\")");
			pw.println("public class "+javaClassName+" extends "+AbstractFilter.class.getName().replace('$', '.')+" {");
			pw.println("  public "+javaClassName+"(final SAMFileHeader header) {");
			pw.println("  super(header);");
			pw.println("  }");
			if(user_code_is_body)
				{
				pw.println("   //user's code starts here");
				pw.println(code);
				pw.println("   // user's code ends here");
				}
			else
				{
				pw.println("  @Override");
				pw.println("  public Object apply(final SAMRecord record) {");
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
			final Constructor<?> ctor=compiledClass.getDeclaredConstructor(SAMFileHeader.class);
			
			
			samFileReader= openSamReader(oneFileOrNull(args));
			final SAMFileHeader header=samFileReader.getFileHeader();
			final AbstractFilter filter = (AbstractFilter)ctor.newInstance(header);
			sw = writingBamArgs.openSAMFileWriter(outputFile,header, true);

			
			long count=0L;
		        final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header).logger(LOG);
		        iter = samFileReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord record=progress.watch(iter.next());
				final Object result = filter.apply(record);
				
				// result is an array of a collection of reads
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
					// write all of reads
					for(final Object item:col)
						{
						if(item==null) throw new JvarkitException.UserError("item in array is null");
						if(!(item instanceof SAMRecord)) throw new JvarkitException.UserError("item in array is not a SAMRecord "+item.getClass());
						++count;
						sw.addAlignment(SAMRecord.class.cast(item));
						if(this.LIMIT>0L && count>=this.LIMIT) break;
						}
					}
				// result is a SAMRecord
				else if(result!=null && (result instanceof SAMRecord)) {
					++count;
					sw.addAlignment(SAMRecord.class.cast(result));
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
						++count;
						failing(record,header);
						}
					else
						{
						sw.addAlignment(record);
						}
					}

				if(this.LIMIT>0L && count>=this.LIMIT) break;
				}
			sw.close();
			/* create empty if never called */
			openFailing(header);
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(samFileReader);
			CloserUtil.close(sw);
			CloserUtil.close(failingReadsWriter);
			}
		}

	public static void main(final String[] args) throws Exception
		{
		new SamJdk().instanceMainWithExit(args);
		}
	}
