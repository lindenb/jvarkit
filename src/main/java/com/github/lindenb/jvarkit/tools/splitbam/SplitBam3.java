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
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;


/**

BEGIN_DOC


Split a BAM by chromosome group. Create EMPTY bams if no reads was found for a given group.
![img](https://chart.googleapis.com/chart?chl=+digraph+G+%7B%0D%0ABWA+-%3E+SPLITBAM%5Blabel%3D%22stdout%22%5D%3B%0D%0A+++SPLITBAM-%3ECHR1_bam%3B%0D%0A+++SPLITBAM-%3ECHR2_bam%3B%0D%0A+++SPLITBAM-%3ECHR3_bam%3B+%0D%0A+++CHR1_bam+-%3E+CHR1_vcf%3B%0D%0A+++CHR2_bam+-%3E+CHR2_vcf%3B%0D%0A+++CHR3_bam+-%3E+CHR3_vcf%3B%0D%0A+++CHR1_vcf+-%3E+merged_vcf%3B%0D%0A+++CHR2_vcf+-%3E+merged_vcf%3B%0D%0A+++CHR3_vcf+-%3E+merged_vcf%3B%0D%0A+%7D%0D%0A++++++++&cht=gv)

file output MUST contain the word '__GROUPID__'



### Example


the content of 'split_g1k_v37_01.txt'



```

CHROMS_01_09	1 2 3 4 5 6 7 8 9
CHROMS_10_0Y	10 11 12 13 14 15 16 17 18 19 20 21 22 X Y 
CHROMS_OTHER	MT GL000207.1 GL000226.1 GL000229.1 GL000231.1 GL000210.1 GL000239.1 GL000235.1 GL000201.1 GL000247.1 GL000245.1 GL000197.1 GL000203.1 GL000246.1 GL000249.1 GL000196.1 GL000248.1 GL000244.1 GL000238.1 GL000202.1 GL000234.1 GL000232.1 GL000206.1 GL000240.1 GL000236.1 GL000241.1 GL000243.1 GL000242.1 GL000230.1 GL000237.1 GL000233.1 GL000204.1 GL000198.1 GL000208.1 GL000191.1 GL000227.1 GL000228.1 GL000214.1 GL000221.1 GL000209.1 GL000218.1 GL000220.1 GL000213.1 GL000211.1 GL000199.1 GL000217.1 GL000216.1 GL000215.1 GL000205.1 GL000219.1 GL000224.1 GL000223.1 GL000195.1 GL000212.1 GL000222.1 GL000200.1 GL000193.1 GL000194.1 GL000225.1 GL000192.1 

```



split the output of bwa sampe on the fly:



```

bwa mem (...) | samtools sort (...) | \
java -jar dist/splitbam3.jar \
	-o TESTSPLITBAM/__GROUPID__.bam \
	-m \
	-g split_g1k_v37_01.txt 


[Fri Jul 26 13:25:56 CEST 2013] Executing as lindenb@master on Linux 2.6.32-358.6.2.el6.x86_64 amd64; OpenJDK 64-Bit Server VM 1.7.0_19-mockbuild_2013_04_17_19_18-b00; Picard version: null
INFO	2013-07-26 13:25:56	SplitBam	reading stdin
INFO	2013-07-26 13:25:56	SplitBam	opening TESTSPLITBAM/CHROMS_01_09.bam
INFO	2013-07-26 13:25:57	SplitBam	opening TESTSPLITBAM/CHROMS_10_0Y.bam
INFO	2013-07-26 13:25:58	SplitBam	opening TESTSPLITBAM/CHROMS_OTHER.bam
INFO	2013-07-26 13:35:58	SplitBam	closing group CHROMS_01_09
INFO	2013-07-26 13:35:59	SplitBam	closing group CHROMS_10_0Y
INFO	2013-07-26 13:35:59	SplitBam	closing group CHROMS_OTHER
INFO	2013-07-26 13:36:00	SplitBam	closing group Unmapped
Runtime.totalMemory()=1916600320

```







END_DOC
*/
@Program(name="splitbam3",description="Split a BAM by chromosome group")
public class SplitBam3 extends Launcher
	{

	private static final Logger LOG = Logger.build(SplitBam3.class).make();


	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-g","--groupfile"},description="Chromosome group file. Interval are 1 based")
	private File chromGroupFile = null;

	@Parameter(names={"-m","--mock"},description="add mock record if no samRecord saved in bam")
	private boolean ADD_MOCK_RECORD = false;

	@Parameter(names={"-u","--unmapped"},description="unmapped chromosome name")
	private String UNDERTERMINED_NAME = "Unmapped";

	//@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
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
					SplitBam3.this.outputFile.getParentFile(),
					SplitBam3.this.outputFile.getName().replaceAll(
							SplitBam3.REPLACE_GROUPID,
							this.groupName
							));
			}
		
		public void open(final SAMFileHeader src)
			{
			SAMFileWriterFactory samFileWriterFactory=writingBamArgs.createSAMFileWriterFactory();
			samFileWriterFactory.setMaxRecordsInRam(writingSortingCollection.getMaxRecordsInRam());
			samFileWriterFactory.setTempDirectory(writingSortingCollection.getTmpDirectories().get(0));

			
			final File fileout=getFile();
			LOG.info("opening BAM file "+fileout);
			final File parent=fileout.getParentFile();
			if(parent!=null) {
				parent.mkdirs();
				samFileWriterFactory.setTempDirectory(writingSortingCollection.getTmpDirectories().get(0));
			}

			
			this.header= src.clone();
			this.header.addComment(
					"Processed with "+getProgramCommandLine()+
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
				throw new IOException("file.is.missing.dict");
				}
			if(!SAMFileHeader.SortOrder.coordinate.equals(srcHeader.getSortOrder()))
				{
				throw new IOException("Input file Bad sort-order: "
						+ srcHeader.getSortOrder()+" exepected coordinate."
						);
				}
			
			this.underminedGroup = new SplitGroup(UNDERTERMINED_NAME);
			this.name2group.put(UNDERTERMINED_NAME,this.underminedGroup);
	
			
			
			if(this.chromGroupFile!=null)
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
	public int doWork(List<String> args) {
		SamReader sfr=null;
		try
			{
			if(outputFile==null || !outputFile.getName().contains(REPLACE_GROUPID))
				{
				LOG.error("output file pattern undefined or doesn't contain "+REPLACE_GROUPID+" : "+this.outputFile);
				return -1;
				}
			if(!outputFile.getName().endsWith(".bam"))
				{
				LOG.error("output file must end with '.bam'");
				return -1;
				}
			
			sfr = super.openSamReader(oneFileOrNull(args));
			
					
			scan(sfr);
			
			LOG.info("Done");
			sfr.close();sfr=null;
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
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
