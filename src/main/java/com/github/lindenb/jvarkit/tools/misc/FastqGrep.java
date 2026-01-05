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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.StringUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.fastq.FastqPairedWriter;
import com.github.lindenb.jvarkit.fastq.FastqRecordPair;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassFastqLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.log.Logger;

/**
BEGIN_DOC

TODO

## Deprecation:

use picard

END_DOC
*/

@Program(
	name="fastqgrep",
	description="Grep reads names in fastq",
	deprecatedMsg="use picard",
	keywords={"fastq"},
	modificationDate="20220208",
	jvarkit_amalgamion = true
	)
public class FastqGrep
	extends OnePassFastqLauncher
	{
	private static final Logger LOG = Logger.of(FastqGrep.class);


	@Parameter(names="-f",description=" file containing a list of read names")
	private Path readNameFile=null;
	@Parameter(names="-R",description="add the read")
	private Set<String> readNamesInput =new HashSet<>();
	@Parameter(names="-n",description="when found, remove the read from the list of names when found more that 'n' time (increase speed)")
	private int n_before_remove=-1;
	@Parameter(names="-V",description="invert)")
	private boolean inverse=false;
	
	private Map<String,Integer> readNames=new HashMap<String,Integer>(); 

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private String getReadName(final FastqRecord r)
		{
		return getReadName(r.getReadName());
		}
	
	private String getReadName(String s)
		{
		final int beg=(s.startsWith(FastqConstants.SEQUENCE_HEADER)?1:0);
		int end=s.indexOf(' ');
		if(end==-1) end=s.length();
		s= s.substring(beg, end);
		return s;
		}
	
	private boolean keepRead(final String readName) {
		boolean keep = false;
		Integer count= this.readNames.get(readName);
		if(count!=null)
			{
			keep=true;
			}
		if(this.inverse) keep=!keep;
		
		
		if(this.n_before_remove!=-1 && !this.inverse && keep)
			{
			count++;
			if(count>=this.n_before_remove)
				{
				this.readNames.remove(readName);
				}
			else
				{
				this.readNames.put(readName,count);
				}
			}
		return keep;
		}
	
	
	@Override
	protected int runSingleEnd(final FastqReader r, FastqWriter out) throws IOException {
		while(r.hasNext() && !this.readNames.isEmpty())
			{
			final FastqRecord fastq=r.next();
			if(keepRead(getReadName(fastq))) {
				out.write(fastq);
				}
			}
		return 0;
		}
	@Override
	protected int runPairedEnd(CloseableIterator<FastqRecordPair> iter, FastqPairedWriter fws) throws IOException {
		while(iter.hasNext() && !this.readNames.isEmpty()) {
			final FastqRecordPair rec = iter.next();
			if(keepRead(getReadName(rec.get(0))) || keepRead(getReadName(rec.get(1)))) {
				fws.write(rec);
				}
			}
		return 0;
		}
	@Override
	protected int beforeFastq() {
		try  {
			if(this.readNameFile!=null)
				{
				try(BufferedReader in=IOUtils.openPathForBufferedReading(this.readNameFile)) {
			    	String line;
			    	while((line=in.readLine())!=null)
			    		{
			    		line=line.trim();
			    		if(StringUtil.isBlank(line)) continue;
			    		this.readNames.put(getReadName(line),0);
			    		}
					}
				}
			
			for(final String r: this.readNamesInput)
				{	
				this.readNames.put(getReadName(r),0);
				}
			
			if(readNames.isEmpty())
	    		{
	    		LOG.warn("no read name found.");
	    		}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new FastqGrep().instanceMainWithExit(args);

	}

}
