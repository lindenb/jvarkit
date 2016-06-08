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
package com.github.lindenb.jvarkit.tools.splitbam;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
//import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.ProgressLoggerInterface;

/***
 * 
 * SplitBam
 *
 */
public class SplitBam3 extends AbstractSplitBam3
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(SplitBam3.class);
	private final static String REPLACE_GROUPID="__GROUPID__";
	private long id_generator=System.currentTimeMillis();
	private java.util.Map<String,SplitGroup> name2group=new java.util.HashMap<String,SplitGroup>();
	private IntervalTreeMap<SplitGroup> interval2group = new IntervalTreeMap<SplitGroup>();
	private SplitGroup underminedGroup=null;
	
	private class SplitGroup
		implements SAMFileWriter
		{
		final String groupName;
		SAMFileHeader header=null;
		SAMFileWriter _writer;
		long count=0L;
		@SuppressWarnings("unused")
		ProgressLoggerInterface progress;
		
		SplitGroup(final String groupName)
			{
			this.groupName=groupName;
			}
		
		@Override
		public SAMFileHeader getFileHeader() {
			return header;
			}
		
		public File getFile()
			{
			return new File(
					SplitBam3.this.getOutputFile().getParentFile(),
					SplitBam3.this.getOutputFile().getName().replaceAll(
							SplitBam3.REPLACE_GROUPID,
							this.groupName
							));
			}
		
		public void open(final SAMFileHeader src)
			{	
			final SAMFileWriterFactory samFileWriterFactory= SplitBam3.this.createSAMFileWriterFactory();
			samFileWriterFactory.setMaxRecordsInRam(maxRecordsInRam);

			
			final File fileout=getFile();
			LOG.info("opening BAM file "+fileout);
			final File parent=fileout.getParentFile();
			if(parent!=null) {
				parent.mkdirs();
				samFileWriterFactory.setTempDirectory(parent);
			}

			
			this.header= src.clone();
			this.header.addComment(
					"Processed with "+getName()+
					" version:"+getVersion()+
					"CommandLine:"+getProgramCommandLine()
					);
			
			this.header.setSortOrder(SortOrder.coordinate);
			samFileWriterFactory.setCreateIndex(true);
				
			
			this._writer = samFileWriterFactory.makeBAMWriter(
				this.header,
				true,
				fileout
				);
			}
		
		@Override
		public void addAlignment(final SAMRecord rec)
			{
			this._writer.addAlignment(rec);
			++this.count;
			}
		
		@Override
		public void setProgressLogger(final ProgressLoggerInterface progress) {
			this.progress=progress;
			}
		
		@Override
		public void close()
			{
			LOG.info("CLOSING "+this.groupName+" N="+this.count);
			if(count==0L && SplitBam3.this.ADD_MOCK_RECORD)
				{
				final List<SAMReadGroupRecord> G=getFileHeader().getReadGroups();
				final String bases="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
				final SAMRecordFactory f=new DefaultSAMRecordFactory();
				++id_generator;
				for(int i=0;i< 2;++i)
					{
					final SAMRecord rec=f.createSAMRecord(getFileHeader());
					rec.setFirstOfPairFlag(i%2==0);
					rec.setSecondOfPairFlag(i%2==1);
					rec.setReadBases(bases.getBytes());
					rec.setMappingQuality(0);
					rec.setBaseQualityString(bases.replace('N', '#'));
					rec.setReadUnmappedFlag(true);
					rec.setMateUnmappedFlag(true);
					rec.setReadPairedFlag(true);
					String readName="MOCKREAD"+(id_generator)+":1:190:289:82";
					rec.setReadName(readName);
					LOG.info("generating mock read: "+readName);
					rec.setAttribute("MK",1);
					if(G!=null && !G.isEmpty())
						{
						rec.setAttribute("RG", G.get(0).getId());
						}
					this.addAlignment(rec);
					}
				}
			CloserUtil.close(this._writer);
			this._writer=null;
			}
		}
	
	private SplitBam3()
		{
		
		}
	
	private SplitGroup getGroupFromInterval(Interval interval)
		{
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

			
		
	
	@SuppressWarnings("resource")
	private void scan(final SamReader samFileReader) throws Exception
		{
		try
			{
			final SAMFileHeader srcHeader= samFileReader.getFileHeader();
			final SAMSequenceDictionary samSequenceDictionary=samFileReader.getFileHeader().getSequenceDictionary();
			if(samSequenceDictionary==null || samSequenceDictionary.isEmpty())
				{
				throw new IOException(getMessageBundle("file.is.missing.dict"));
				}
			if(!SAMFileHeader.SortOrder.coordinate.equals(srcHeader.getSortOrder()))
				{
				throw new IOException(getMessageBundle("Input file Bad sort-order: "
						+ srcHeader.getSortOrder()+" exepected coordinate."
						));
				}
			
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
							SAMSequenceRecord ssr=samSequenceDictionary.getSequence(segment);
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
			
		
			/* open all output bams */
			for(final SplitGroup g:this.name2group.values())
				{
				g.open(srcHeader);
				}
			
	        
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samFileReader.getFileHeader()==null?null:samFileReader.getFileHeader().getSequenceDictionary());
	        
			for(Iterator<SAMRecord> iter=samFileReader.iterator();
					iter.hasNext(); )
				{
				final SAMRecord record = progress.watch(iter.next());
			
				Interval interval=null;
				if( record.getReadUnmappedFlag() )
					{
					if(record.getReadPairedFlag() && !record.getMateUnmappedFlag())
						{
						interval= new Interval(
								record.getMateReferenceName(),
								record.getMateAlignmentStart(),
								record.getMateAlignmentStart()
								);
						}
					}
				else
					{
					interval= new Interval(
							record.getReferenceName(),
							record.getAlignmentStart(),
							record.getAlignmentStart()
							);
					
					}
				SplitGroup splitGroup = 
						(
						interval == null ?
						null :
						this.getGroupFromInterval(interval)
						);
				
				if(splitGroup==null) splitGroup=this.underminedGroup;
				
				splitGroup.addAlignment(record);
				}
			
			samFileReader.close();
			
			/* copenlose all */
			for(final SplitGroup g:this.name2group.values())
				{
				g.close();
				}
			
			progress.finish();
			}
		catch(final Exception error)
			{
			LOG.error("failure:",error);
			for(final SplitGroup g:this.name2group.values())
				{
				g.close();
				final File f=g.getFile();
				if(f.exists())
					{
					LOG.info("Delete "+f);
					f.delete();
					}
				}
			throw error;
			}
		}
	
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		SamReader sfr=null;
		try
			{
			if(getOutputFile()==null || !getOutputFile().getName().contains(REPLACE_GROUPID))
				{
				return wrapException("output file pattern undefined or doesn't contain "+REPLACE_GROUPID+" : "+this.getOutputFile());
				}
			if(!getOutputFile().getName().endsWith(".bam"))
				{
				return wrapException("output file must end with '.bam'");
				}
			
			sfr = super.openSamReader(inputName);
			
					
			scan(sfr);
			
			LOG.info("Done");
			sfr.close();sfr=null;
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(sfr);
			}
		
		}
	
	
	
	public static void main(String[] args) throws Exception
		{
		new SplitBam3().instanceMainWithExit(args);
		}
	
	}
