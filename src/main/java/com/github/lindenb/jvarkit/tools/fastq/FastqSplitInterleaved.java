/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.fastq;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

/**
 * FastqSplitInterleaved
 * @author lindenb
 *
 */
public class FastqSplitInterleaved extends AbstractCommandLineProgram
	{
	private FastqSplitInterleaved()
		{
		
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FastqSplitInterleaved";
		}
	
	@Override
	public String getProgramDescription() {
		return "Split interleaved Fastq files.";
		}
	
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -a (fastq1 file or '-' for stdout). Ignore 1st read if omitted. Optional.");
		out.println(" -b (fastq2 file or '-' for stdout). Ignore 2nd read if omitted. Optional.");
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{
		String fileout[]={null,null};
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"a:b:"))!=-1)
			{
			switch(c)
				{
				case 'a': fileout[0]= opt.getOptArg();break;
				case 'b': fileout[1]= opt.getOptArg();break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		FastqReader r1=null;
		FastqWriter writers[]={null,null};
		try
			{
			if(opt.getOptInd() == args.length)
				{
				r1=new FourLinesFastqReader(System.in);				
				}
			else if(opt.getOptInd()+1==args.length)
				{
				r1=new FourLinesFastqReader(new File(args[opt.getOptInd()]));				
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			
			if(fileout[0]==null && fileout[1]==null)
				{
				error("Both outputs are undefined.");
				return -1;
				}
			
			for(int i=0;i<2;++i)
				{
				if(fileout[i]==null)
					{
					writers[i]=new BasicFastqWriter(new PrintStream(new NullOuputStream()));
					}
				else if(fileout[i].equals("-"))
					{
					if(i==1 && "-".equals(fileout[0]))
						{
						writers[i]=writers[0];
						}
					else
						{
						writers[i]=new BasicFastqWriter(System.out);
						}
					}
				else
					{
					if(i==1 && fileout[1].equals(fileout[0]))
						{
						writers[i]=writers[0];
						}
					else
						{
						writers[i]=new BasicFastqWriter(new File(fileout[i]));
						}
					}
				}

			
			
			FastqRecord records[]={null,null};
			
			while(r1.hasNext())
				{
				records[0]=r1.next();
				if(!r1.hasNext())
					{
					throw new IOException(getMessageBundle("fastq.paired.read.missing"));
					}
				records[1]=r1.next();
				
				for(int i=0;i< 2;++i)
					{
					writers[i].write(records[i]);
					}
				
				}
			if(r1.hasNext())
				{
				throw new IOException("Illegal number of reads in fastq");
				}
			
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r1);
			CloserUtil.close(writers[0]);
			CloserUtil.close(writers[1]);
			}
		}
	
	public static void main(String[] args)
		{
		new FastqSplitInterleaved().instanceMainWithExit(args);
		}

}
