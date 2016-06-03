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
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.Collection;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class SplitVcf
	extends AbstractSplitVcf
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SplitVcf.class);

	private final static String REPLACE_GROUPID="__GROUPID__";
	private java.util.Map<String,SplitGroup> name2group=new java.util.HashMap<String,SplitGroup>();
	private IntervalTreeMap<SplitGroup> interval2group = new IntervalTreeMap<SplitGroup>();
	private SplitGroup underminedGroup=null;

	
	private class SplitGroup implements Closeable
		{
		final String groupName;
		VCFHeader header=null;
		VariantContextWriter _writer;
		
		SplitGroup(final String groupName)
			{
			this.groupName=groupName;
			}
		

		@Override
		public void close() {
			CloserUtil.close(_writer);
			
		}

		
		public File getFile()
			{
			return new File(
					SplitVcf.this.OUT_FILE_PATTERN.replaceAll(
							SplitVcf.REPLACE_GROUPID,
							this.groupName
							));
			}
		
		public void open(final VCFHeader src)
			{	
	
			
			final File fileout=getFile();
			LOG.info("opening BAM file "+fileout);
			final File parent=fileout.getParentFile();
			if(parent!=null) {
				parent.mkdirs();
			}
	
			
			this.header= new VCFHeader(src);
			try {
				this._writer = VCFUtils.createVariantContextWriter(fileout);
			} catch (IOException e) {
				throw new RuntimeIOException(e);
			}
			this._writer.writeHeader(this.header);
			}
		}
		
		public SplitVcf()
			{
			}
		
		
		
		private SplitGroup getGroupFromInterval(Interval interval) {
			final Collection<SplitGroup> groups=this.interval2group.getOverlapping(interval);
			if(groups==null || groups.isEmpty()) return null;
			if(groups.size()!=1) throw new IllegalStateException();
			return groups.iterator().next();
			}
	
	private void put(final String groupName,final Interval interval)
		{
		SplitGroup splitgroup = this.getGroupFromInterval(interval);
		if(splitgroup!=null && !splitgroup.groupName.equals(groupName))
			{
			throw new IllegalArgumentException("chrom "+interval+" already used in "+splitgroup.groupName);
			}
		splitgroup = name2group.get(groupName);
				
		if(splitgroup==null)
			{
			splitgroup = new SplitGroup(groupName);
			
			this.name2group.put(groupName,splitgroup);
			}
		this.interval2group.put(interval, splitgroup);
		}
	
		
			 
			
			@Override
			protected Collection<Throwable> call(String inputName) throws Exception {
				if(!OUT_FILE_PATTERN.contains(REPLACE_GROUPID))
					{
					return wrapException("output file pattern undefined or doesn't contain "+REPLACE_GROUPID+" : "+this.OUT_FILE_PATTERN);
					}
			if(!(OUT_FILE_PATTERN.endsWith(".vcf") || OUT_FILE_PATTERN.endsWith(".vcf.gz")))
				{
				return wrapException("output file must end with '.vcf' or '.vcf.gz'");
				}
			
			VcfIterator in=openVcfIterator(inputName);
			final SAMSequenceDictionary samSequenceDictionary=in.getHeader().getSequenceDictionary();
			
			this.underminedGroup = new SplitGroup(UNDERTERMINED_NAME);
			this.name2group.put(UNDERTERMINED_NAME,this.underminedGroup);
	
			
			
			if(super.chromGroupFile!=null)
				{
				BufferedReader r=IOUtils.openFileForBufferedReading(chromGroupFile);
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					final String tokens[] =line.split("[ \t,]+");
					final String groupName=tokens[0].trim();
					if(groupName.isEmpty()) throw new IOException("Empty group name in "+line);
					if(this.UNDERTERMINED_NAME.equals(groupName))  throw new IOException("Group cannot be named "+UNDERTERMINED_NAME);
					if(this.name2group.containsKey(groupName))  throw new IOException("Group defined twice "+groupName);
					for(int i=1;i< tokens.length;i++)
						{
						String sequence;
						int start;
						int end;
						String segment = tokens[i].trim();
						
						if(segment.isEmpty()) continue;
						
						int colon= segment.indexOf(':');
						if(colon==-1)
							{
							final SAMSequenceRecord ssr=samSequenceDictionary.getSequence(segment);
							if(ssr==null)
								{
								throw new IOException("Unknown chromosome , not in dict \""+segment+"\"");
								}
							sequence = segment;
							start = 1;
							end = ssr.getSequenceLength();
							}
						else
							{
							int hyphen  = segment.indexOf('-',colon);
							if(hyphen==-1)  throw new IOException("Bad segment:"+segment);
							sequence = segment.substring(0,colon);
							if(samSequenceDictionary.getSequence(sequence)==null)
								 throw new IOException("Unknown chromosome , not in dict "+
										 segment);
							
							//+1 because interval are 1-based
							start = 1+Integer.parseInt(segment.substring(colon+1,hyphen));
							end = Integer.parseInt(segment.substring(hyphen+1));
							}
						
						final Interval interval = new Interval(sequence, start, end);
						this.put(groupName,interval);
						}
					}
				r.close();
				}
			else
				{
				LOG.info("creating default split interval");
				for(final SAMSequenceRecord seq:samSequenceDictionary.getSequences())
					{
					final String groupName=seq.getSequenceName();
					final Interval interval= new Interval(groupName, 1, seq.getSequenceLength());
					this.put(groupName,interval);
					}
				}
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samSequenceDictionary);
	        
			while(in.hasNext())
				{
				final VariantContext record = progress.watch(in.next());
			
				Interval interval=new Interval(
						record.getContig(),
						record.getStart(),
						record.getEnd()
						);
				SplitGroup splitGroup = this.getGroupFromInterval(interval);
				
				if(splitGroup==null) splitGroup=this.underminedGroup;
				
				splitGroup._writer.add(record);
				}
			
			progress.finish();
			in.close();

		
			/* open all output vcf */
			for(final SplitGroup g:this.name2group.values())
				{
				g.open(in.getHeader());
				}
			for(final SplitGroup g:this.name2group.values())
				{
				g.close();
				}
			return RETURN_OK;
			}
		 	
		
		public static void main(String[] args)
			{
			new SplitVcf().instanceMainWithExit(args);
			}

		}
