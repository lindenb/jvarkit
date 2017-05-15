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
	
	
	private PadEmptyFastq()
		{
		}
	
	
	private void copyTo(FastqReader r,FastqWriter w)
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
			if(rec.getReadString().isEmpty())
				{
				++nFill;
				if(padLength<1)
					{
					padLength=DEFAULT_LENGTH;
					}
				if(fillN==null)
					{
					StringBuilder b1=new StringBuilder();
					while(b1.length()< padLength) b1.append("N");
					fillN=b1.toString();
					fillQ=fillN.replace('N', '#');
					}
				
				rec=new FastqRecord(
						rec.getReadHeader(),
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
	public int doWork(List<String> args) {
		
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
				FastqReader fqr=new FourLinesFastqReader(System.in);
				copyTo(fqr,fqw);
				fqr.close();
				}
			else
				{
				for(String filename:args)
					{
					LOG.info("Reading from "+filename);
					FastqReader fqr=new FourLinesFastqReader(new File(filename));
					copyTo(fqr,fqw);
					fqr.close();
					}
				}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(fqw);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new PadEmptyFastq().instanceMainWithExit(args);
		}

}
