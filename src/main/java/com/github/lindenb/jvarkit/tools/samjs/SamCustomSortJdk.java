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
package com.github.lindenb.jvarkit.tools.samjs;

import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.lang.reflect.Constructor;
import java.util.Comparator;
import java.util.Date;
import java.util.List;
import java.util.Random;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.OpenJdkCompiler;
import com.github.lindenb.jvarkit.tools.misc.IlluminaReadName;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iso8601Date;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;


/**
BEGIN_DOC

## Motivation

Sort a BAM using a java expression compiled at runtime

## How it works

The user provides some java code of a `ComparatorSamRecord>`. This class extends `SamCustomSortJdk.AbstractSamComparator`.
At the time of writting, the  `SamCustomSortJdk.AbstractSamComparator` is:


```java
public static abstract class AbstractSamComparator
	implements Comparator<SAMRecord>
	{
	//input SAM header
	protected final SAMFileHeader header;
	//coordinate comparator 
	private final SAMRecordCoordinateComparator coordinateComparator;
	// query name comparator
	private final SAMRecordQueryNameComparator queryNameComparator;
	protected AbstractSamComparator(final SAMFileHeader header) {
		this.header = header;
		this.coordinateComparator = new SAMRecordCoordinateComparator();
		this.queryNameComparator = new SAMRecordQueryNameComparator();
		}
	// get a SAMRecordCoordinateComparator ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordCoordinateComparator.html ), never null 
	protected SAMRecordCoordinateComparator getCoordinateComparator() {
		return this.coordinateComparator;
	}
	// get a SAMRecordQueryNameComparator, never null 
	protected SAMRecordQueryNameComparator getQueryNameComparator() {
		return this.queryNameComparator;
	}
	// sort reads using 'SAMRecordCoordinateComparator.fileOrderCompare' https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordCoordinateComparator.html
	protected int fileOrderCompare(final SAMRecord R1, final SAMRecord R2) {
    	return getCoordinateComparator().fileOrderCompare(R1, R2);
    	}
	}
```
 
The user's code can be the code of the function `int compareTo(final SAMRecord R1,final SAMRecord R2)`, for example to sort on the mapping quality:

```java
return R1.getMappingQuality() - R2.getMappingQuality();
```

or the whole body of the comparator (option --body )


```java
private int mapq(final SAMRecord R) {
return R.getMappingQuality();
}

public int compareTo(final SAMRecord R1,final SAMRecord R2) {
return mapq(R1) - mapq(R2);
}
```

## Example


sort on amount of reference sequence covered, using the cigar string

```
 java -jar dist/samcustomsortjdk.jar  --body -e 'private int score(final SAMRecord R) { if(R.getReadUnmappedFlag() || R.getCigar()==null) return 0; return R.getCigar().getReferenceLength();} @Override public int compare(SAMRecord a,SAMRecord b) { return Integer.compare(score(a),score(b));}' in.bam
 ```

## History

 * 2019-01: migrating to openjkd: switched to in-memory compiling to external compiling

END_DOC
*/
@Program(name="samcustomsortjdk",
	description="Sort a BAM file using a java expression compiled at runtime.",
	keywords={"sam","bam","java","jdk","sort"},
	biostars={305181}
	)
public class SamCustomSortJdk
	extends Launcher
	{
	private static final Logger LOG = Logger.build(SamCustomSortJdk.class).make();
	@SuppressWarnings("unused")
	private IlluminaReadName _fool_javac = null;
	
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	@Parameter(names={"-e","--expression"},description="java expression")
	private String scriptExpr=null;
	@Parameter(names={"-f","--file"},description="java file. Either option -e or -f is required.")
	private File scriptFile =null;	
	@Parameter(names={"--nocode"},description=" Don't show the generated code")
	private boolean hideGeneratedCode=false;
	@Parameter(names={"--body"},description="user's code is the whole body of the filter class, not just the 'apply' method.")
	private boolean user_code_is_body=false;
	@Parameter(names={"--saveCodeInDir"},description="Save the generated java code in the following directory")
	private File saveCodeInDir=null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	public static abstract class AbstractSamComparator
		implements Comparator<SAMRecord>
		{
		/** input SAM header */
		protected final SAMFileHeader header;
		/** coordinate comparator */
		private final SAMRecordCoordinateComparator coordinateComparator;
		/** query name comparator */
		private final SAMRecordQueryNameComparator queryNameComparator;
		
		protected AbstractSamComparator(final SAMFileHeader header) {
			this.header = header;
			this.coordinateComparator = new SAMRecordCoordinateComparator();
			this.queryNameComparator = new SAMRecordQueryNameComparator();
			}
		/** get a SAMRecordCoordinateComparator ( https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordCoordinateComparator.html ), never null */
		protected SAMRecordCoordinateComparator getCoordinateComparator() {
			return this.coordinateComparator;
		}
		/** get a SAMRecordQueryNameComparator, never null */
		protected SAMRecordQueryNameComparator getQueryNameComparator() {
			return this.queryNameComparator;
		}
		/** sort reads using 'SAMRecordCoordinateComparator.fileOrderCompare' https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/SAMRecordCoordinateComparator.html */
		protected int fileOrderCompare(final SAMRecord R1, final SAMRecord R2) {
	    	return getCoordinateComparator().fileOrderCompare(R1, R2);
	    	}
		
		}

	
	private static class StableSort implements SAMRecordComparator
		{
		private final SAMRecordCoordinateComparator coordinateComparator;
		final  Comparator<SAMRecord> delegate;
		StableSort(final  Comparator<SAMRecord> delegate) {
			this.delegate = delegate;
			this.coordinateComparator = new SAMRecordCoordinateComparator();
		}
		@Override
		public int compare(final SAMRecord o1,final SAMRecord o2) {
			final int i= delegate.compare(o1, o2);
			if(i!=0) return i;
			return fileOrderCompare(o1, o2);
			}
		@Override
		public int fileOrderCompare(final SAMRecord R1, final SAMRecord R2) {
	    	return coordinateComparator.fileOrderCompare(R1, R2);
	    	}
		}

	
	
	public SamCustomSortJdk()
		{
		}

	
	@Override
	public int doWork(final List<String> args) {
		SAMRecordIterator iter=null;
		SamReader samFileReader=null;
		SAMFileWriter sw=null;
		SortingCollection<SAMRecord> sorter=null;
		CloseableIterator<SAMRecord> iter2=null;
		try
			{
			final String code;
			
			if(this.scriptFile!=null)
				{
				code = IOUtil.slurp(this.scriptFile);
				}
			else if(!StringUtil.isBlank(this.scriptExpr))
				{
				code = this.scriptExpr;
				}
			else
				{
				LOG.error("Option -e or -f are required. The content of those empty mut be not empty");
				return -1;
				}

			final Random rand= new  Random(System.currentTimeMillis());
			final String javaClassName =SamCustomSortJdk.class.getSimpleName()+
					"Custom"+ Math.abs(rand.nextInt());
			
			final StringWriter codeWriter=new StringWriter();
			final PrintWriter pw = new PrintWriter(codeWriter);
			pw.println("import java.util.*;");
			pw.println("import java.util.stream.*;");
			pw.println("import java.util.function.*;");
			pw.println("import htsjdk.samtools.*;");
			pw.println("import htsjdk.samtools.util.*;");
			pw.println("import com.github.lindenb.jvarkit.tools.misc.IlluminaReadName;");
			pw.println("import javax.annotation.processing.Generated;");

			pw.println("@Generated(value=\""+SamCustomSortJdk.class.getSimpleName()+"\",date=\""+ new Iso8601Date(new Date()) +"\")");
			pw.println("public class "+javaClassName+" extends "+
					AbstractSamComparator.class.getName().replace('$', '.')+" {");
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
				pw.println("  public int compare(final  SAMRecord R1,final SAMRecord R2) {");
				pw.println("   /** user's code starts here */");
				pw.println(code);
				pw.println(    "/** user's code ends here */");
				pw.println("   }");
				}
			pw.println("}");
			pw.flush();
			
			
			if(!hideGeneratedCode)
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
					LOG.error(err);
					return -1;
					}
				finally
					{
					CloserUtil.close(cw);
					}
				}
			
			final OpenJdkCompiler inMemoryCompiler = OpenJdkCompiler.getInstance();
			final Class<?> compiledClass = inMemoryCompiler.compileClass(
					javaClassName,
					codeWriter.toString()
					);
			final Constructor<?> ctor=compiledClass.getDeclaredConstructor(SAMFileHeader.class);

			
			samFileReader= openSamReader(oneFileOrNull(args));
			final SAMFileHeader headerIn = samFileReader.getFileHeader();
			@SuppressWarnings("unchecked")
			final StableSort customComparator = new StableSort(( Comparator<SAMRecord>)ctor.newInstance(headerIn));
			final BAMRecordCodec bamRecordCodec=new BAMRecordCodec(headerIn);
			
			sorter =SortingCollection.newInstance(
						SAMRecord.class,
						bamRecordCodec,
						customComparator,
						this.writingSortingCollection.getMaxRecordsInRam(),
						this.writingSortingCollection.getTmpPaths()
						);
			sorter.setDestructiveIteration(true);
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(headerIn).logger(LOG);
			iter = samFileReader.iterator();
			while(iter.hasNext())
				{
				sorter.add( progress.watch(iter.next()));
				}
			samFileReader.close();samFileReader=null;
			sorter.doneAdding();
			
			
			
			
			final SAMFileHeader headerOut=headerIn.clone();
			headerOut.setSortOrder(SAMFileHeader.SortOrder.unsorted);
			headerOut.addComment(getProgramName()+" "+getVersion()+" "+getProgramCommandLine());
	        sw = this.writingBamArgs.openSAMFileWriter(this.outputFile,headerOut, false);

	        progress=new SAMSequenceDictionaryProgress(headerIn).logger(LOG);
			iter2 = sorter.iterator();
			while(iter2.hasNext())
				{
				sw.addAlignment( progress.watch(iter2.next()));
				}
			iter2.close();iter2=null;
			sw.close();sw=null;
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
			try { if(sorter!=null) sorter.cleanup();} catch(Exception e) {}
			CloserUtil.close(iter);
			CloserUtil.close(iter2);
			CloserUtil.close(samFileReader);
			CloserUtil.close(sw);
			}
		}

	public static void main(final String[] args) throws Exception
		{
		new SamCustomSortJdk().instanceMainWithExit(args);
		}
	}
