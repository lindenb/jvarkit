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
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
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
public class SplitBam3 extends AbstractCommandLineProgram
	{
	private int maxRecordsInRam=1000000;
	private boolean ADD_MOCK_RECORD=false;
	private final static String REPLACE_GROUPID="__GROUPID__";
	private String OUT_FILE_PATTERN="";
	private String UNDERTERMINED_NAME="Unmapped";
	
	private long id_generator=System.currentTimeMillis();
	private java.util.Map<String,SplitGroup> name2group=new java.util.HashMap<String,SplitGroup>();
	private IntervalTreeMap<SplitGroup> interval2group = new IntervalTreeMap<SplitGroup>();
	private SplitGroup underminedGroup=null;
	
	private class SplitGroup
		implements SAMFileWriter
		{
		String groupName;
		SAMFileHeader header=null;
		SAMFileWriter _writer;
		long count=0L;
		@SuppressWarnings("unused")
		ProgressLoggerInterface progress;
		
		SplitGroup(String groupName)
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
					SplitBam3.this.OUT_FILE_PATTERN.replaceAll(
							SplitBam3.REPLACE_GROUPID,
							this.groupName
							));
			}
		
		public void open(SAMFileHeader src)
			{	
			SAMFileWriterFactory samFileWriterFactory=new SAMFileWriterFactory();
			
			samFileWriterFactory.setTempDirectory(  SplitBam3.this.getTmpDirectories().get(0) );
			samFileWriterFactory.setMaxRecordsInRam(maxRecordsInRam);

			
			File fileout=getFile();
			info("opening BAM file "+fileout);
			File parent=fileout.getParentFile();
			if(parent!=null) parent.mkdirs();

			
			this.header=src.clone();
			this.header.addComment(
					"Processed with "+getProgramName()+
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
		public void addAlignment(SAMRecord rec)
			{
			this._writer.addAlignment(rec);
			++this.count;
			}
		
		@Override
		public void setProgressLogger(ProgressLoggerInterface progress) {
			this.progress=progress;
			}
		
		@Override
		public void close()
			{
			info("CLOSING "+this.groupName+" N="+this.count);
			if(count==0L && SplitBam3.this.ADD_MOCK_RECORD)
				{
				List<SAMReadGroupRecord> G=getFileHeader().getReadGroups();
				String bases="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
				SAMRecordFactory f=new DefaultSAMRecordFactory();
				++id_generator;
				for(int i=0;i< 2;++i)
					{
					SAMRecord rec=f.createSAMRecord(getFileHeader());
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
					info("generating mock read: "+readName);
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
		Collection<SplitGroup> groups=this.interval2group.getOverlapping(interval);
		if(groups==null || groups.isEmpty()) return null;
		if(groups.size()!=1) throw new IllegalStateException();
		return groups.iterator().next();
		}
	
	private void put(String groupName,Interval interval)
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
	private void scan(SamReader samFileReader,File chromGroupFile) throws Exception
		{
		try
			{
			final SAMFileHeader srcHeader= samFileReader.getFileHeader();
			SAMSequenceDictionary samSequenceDictionary=samFileReader.getFileHeader().getSequenceDictionary();
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
	
			
			
			if(chromGroupFile!=null)
				{
				
				
				BufferedReader r=IOUtils.openFileForBufferedReading(chromGroupFile);
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					String tokens[] =line.split("[ \t,]+");
					String groupName=tokens[0].trim();
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
								throw new IOException("Unknown chromosome , not in dict");
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
						
						Interval interval=new Interval(sequence, start, end);
						this.put(groupName,interval);
						}
					}
				r.close();
				}
			else
				{
				for(SAMSequenceRecord seq:samSequenceDictionary.getSequences())
					{
					String groupName=seq.getSequenceName();
					Interval interval= new Interval(groupName, 1, seq.getSequenceLength());
					this.put(groupName,interval);
					}
				}
			
		
			/* open all output bams */
			for(SplitGroup g:this.name2group.values())
				{
				g.open(srcHeader);
				}
			
	        
	       SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samFileReader.getFileHeader()==null?null:samFileReader.getFileHeader().getSequenceDictionary());
	        
			for(Iterator<SAMRecord> iter=samFileReader.iterator();
					iter.hasNext(); )
				{
				SAMRecord record = progress.watch(iter.next());
			
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
			for(SplitGroup g:this.name2group.values())
				{
				g.close();
				}
			
			progress.finish();
			}
		catch(Exception error)
			{
			for(SplitGroup g:this.name2group.values())
				{
				g.close();
				File f=g.getFile();
				if(f.exists())
					{
					info("Delete "+f);
					f.delete();
					}
				}
			throw error;
			}
		}
	
	@Override
	public String getProgramDescription() {
		return "Split a BAM by chromosome group.";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return DEFAULT_WIKI_PREFIX+"SplitBam";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-p (file) output file pattern. MUST contain "+REPLACE_GROUPID+" and end with .bam");
		out.println("-g (file) Chromosome group file. Optional");
		out.println("-u (name) unmapped chromosome name. Optional. Default "+UNDERTERMINED_NAME);
		out.println("-m add mock record if no samRecord saved in bam");
		out.println("-R (int)  max records in RAM "+getMessageBundle("max.records.in.ram"));
		super.printOptions(out);
		}

	@Override
	public int doWork(String[] args)
		{

		File chromGroupFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"u:g:p:mR:"))!=-1)
			{
			switch(c)
				{	
				case 'R': maxRecordsInRam= Integer.parseInt(opt.getOptArg());break;
				case 'u': UNDERTERMINED_NAME= opt.getOptArg();break;
				case 'p': OUT_FILE_PATTERN= opt.getOptArg();break;
				case 'g': chromGroupFile=new File(opt.getOptArg());break;
				case 'm': ADD_MOCK_RECORD=true;break;
				default:
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		SamReader sfr=null;
		try
			{
			
			
			if(!OUT_FILE_PATTERN.contains(REPLACE_GROUPID))
				{
				error("output file pattern undefined or doesn't contain "+REPLACE_GROUPID+" : "+this.OUT_FILE_PATTERN);
				return -1;
				}
			if(!OUT_FILE_PATTERN.endsWith(".bam"))
				{
				error("output file must end with '.bam'");
				return -1;
				}
			
			
			SamReaderFactory samReaderFactory=SamReaderFactory.makeDefault().
					validationStringency(ValidationStringency.SILENT)
					;
			SamInputResource samInputResource=null;	
			 
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				samInputResource = SamInputResource.of(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				samInputResource = SamInputResource.of(new File(filename));
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			sfr=samReaderFactory.open(samInputResource);
					
			scan(sfr,chromGroupFile);
			
			info("Done");
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
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
