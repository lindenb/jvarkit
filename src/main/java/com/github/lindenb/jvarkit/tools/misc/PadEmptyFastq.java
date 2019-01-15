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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.List;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;
/*
BEGIN_DOC

### AWK

compiled version of awk:

```awk
{if(NR%4==2 && length($0)==0) { printf("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n");} else if(NR%4==0 && length($0)==0) { printf("##################################################\n");} else {print;}}
```


### Example

```bash
$ cutadapt -a AGATCGGAAGAGCGTCGT  2> /dev/null in.fastq.gz |\
  java -jar dist/pademptyfastq.jar -o pad.fastq.gz
```

END_DOC
 
 */

@Program(name="pademptyfastq",
description= "Pad empty fastq sequence/qual with N/#",
deprecatedMsg="use awk",
keywords="fastq"
		)
public class PadEmptyFastq extends Launcher
	{
	private static final Logger LOG=Logger.build(PadEmptyFastq.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
    private File outFile=null;
	private static final int DEFAULT_LENGTH=50;
	@Parameter(names={"-N"},description="number of bases/qual to be added.  -1=length of the first read ")
	private int N=-1;//default , will use the first read length
	
	public PadEmptyFastq()
		{
		}
	
	
	private void copyTo(final FastqReader r,final FastqWriter w)
		{
		int padLength=this.N;
		long nReads=0L;
		long nFill=0L;
		String fillN=null;
		String fillQ=null;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			FastqRecord rec=r.next();
			
			if(++nReads%1E6==0)
				{
				LOG.info("Read "+nReads +" reads. empty reads="+nFill);
				}
			if(StringUtils.isBlank(rec.getReadString()))
				{
				++nFill;
				if(padLength<1)
					{
					padLength=DEFAULT_LENGTH;
					}
				if(fillN==null)
					{
					fillN = StringUtils.repeat(padLength, 'N');
					fillQ = StringUtils.repeat(padLength, '#');
					}
				
				rec=new FastqRecord(
						rec.getReadName(),
						fillN,
						rec.getBaseQualityHeader(),
						fillQ
						);
				}
			else if(padLength<1)
				{
				padLength=rec.getReadString().length();
				}
			w.write(rec);
			}
		LOG.info("Done. Read "+nReads +" reads. empty reads="+nFill);
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		FastqWriter fqw=null;		
		try
			{
			if(this.outFile==null)
				{
				LOG.info("writing to stdout");
				fqw=new BasicFastqWriter(stdout());
				}
			else
				{
				LOG.info("writing to "+this.outFile);
				fqw=new FastqWriterFactory().newWriter(this.outFile);
				}
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				final FastqReader fqr=new FourLinesFastqReader(stdin());
				copyTo(fqr,fqw);
				fqr.close();
				}
			else
				{
				for(final String filename:args)
					{
					LOG.info("Reading from "+filename);
					final FastqReader fqr=new FourLinesFastqReader(new File(filename));
					copyTo(fqr,fqw);
					fqr.close();
					}
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(fqw);
			}
		}
	
	public static void main(final String[] args)
		{
		new PadEmptyFastq().instanceMainWithExit(args);
		}

}
