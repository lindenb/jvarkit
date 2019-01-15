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

import java.io.BufferedReader;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.fastq.BasicFastqWriter;
import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;

@Program(
	name="fastqgrep",
	description="Grep reads names in fastq",
	deprecatedMsg="use picard",
	keywords={"fastq"}
	)
public class FastqGrep
	extends Launcher
	{
	private static final Logger LOG = Logger.build(FastqGrep.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names="-f",description=" file containing a list of read names")
	private File readNameFile=null;
	@Parameter(names="-R",description="add the read")
	private Set<String> readNamesInput =new HashSet<>();
	@Parameter(names="-n",description="when found, remove the read from the list of names when found more that 'n' time (increase speed)")
	private int n_before_remove=-1;
	@Parameter(names="-V",description="invert)")
	private boolean inverse=false;
	
	private Map<String,Integer> readNames=new HashMap<String,Integer>(); 

	
	
	
	private FastqGrep()
		{
		}
	
	

	
	private String getReadName(final FastqRecord r)
		{
		return getReadName(r.getReadName());
		}
	
	private String getReadName(String s)
		{
		int beg=(s.startsWith(FastqConstants.SEQUENCE_HEADER)?1:0);
		int end=s.indexOf(' ');
		if(end==-1) end=s.length();
		s= s.substring(beg, end);
		return s;
		}
	private void run(final FastqReader r,final FastqWriter out)
		{
		long nRec=0L;
		r.setValidationStringency(ValidationStringency.LENIENT);
		while(r.hasNext())
			{
			FastqRecord fastq=r.next();
			boolean keep=false;
			String readName=getReadName(fastq);
			Integer count=readNames.get(readName);
			if(count!=null)
				{
				keep=true;
				}
			if(inverse) keep=!keep;
			if(keep)
				{
				++nRec;
				out.write(fastq);
				}
			
			if(n_before_remove!=-1 && !inverse && keep)
				{
				count++;
				if(count>=n_before_remove)
					{
					readNames.remove(readName);
					if(readNames.isEmpty()) break;
					}
				else
					{
					readNames.put(readName,count);
					}
				}
			
			
			}
		LOG.info("Done. N-Reads:"+nRec);
		}
	
	@Override
	public int doWork(List<String> args) {
		BufferedReader in=null;
		FastqWriter out=null;
		try 
			{
			if(this.readNameFile!=null)
				{
				in=IOUtils.openFileForBufferedReading(this.readNameFile);
		    	String line;
		    	while((line=in.readLine())!=null)
		    		{
		    		line=line.trim();
		    		if(line.isEmpty()) continue;
		    		this.readNames.put(getReadName(line),0);
		    		}
		    	in.close();
				}
			
			for(final String r: this.readNamesInput)
				{	
				this.readNames.put(getReadName(r),0);
				}
			
			if(readNames.isEmpty())
	    		{
	    		LOG.warn("no read name found.");
	    		}

			
			if(this.outputFile!=null)
				{
				LOG.info("Writing to "+this.outputFile);
				out=new BasicFastqWriter(this.outputFile);
				}
			else
				{
				LOG.info("Writing to stdout");
				out=new BasicFastqWriter(stdout());
				}
			
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				FastqReader fqR=new FourLinesFastqReader(System.in);
				run(fqR,out);
				fqR.close();
				}
			else for(String fname:args)
				{
				File f=new File(fname);
				LOG.info("Reading from "+f);
				FastqReader fqR=new FourLinesFastqReader(f);
				run(fqR,out);
				fqR.close();
				}
			CloserUtil.close(out);	
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}
	
	public static void main(String[] args) {
		new FastqGrep().instanceMainWithExit(args);

	}

}
