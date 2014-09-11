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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.PicardException;

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

/***
 * 
 * SplitBam
 *
 */
public class SplitBam2 extends AbstractCommandLineProgram
	{
	private boolean GENERATE_EMPTY_BAM=false;
	private boolean ADD_MOCK_RECORD=false;
	private final static String REPLACE_CHROM="__CHROM__";
	private String OUT_FILE_PATTERN="";
	private String UNDERTERMINED_NAME="Unmapped";
	
	private SAMFileWriterFactory samFileWriterFactory=new SAMFileWriterFactory();
	private long id_generator=System.currentTimeMillis();
	

	
	private SplitBam2()
		{
		
		}
	
	private void addMockPair(
			SAMFileWriter sw,
			SAMFileHeader header
			) throws IOException
		{
		List<SAMReadGroupRecord> G=header.getReadGroups();
		String bases="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
		SAMRecordFactory f=new DefaultSAMRecordFactory();
		++id_generator;
		for(int i=0;i< 2;++i)
			{
			SAMRecord rec=f.createSAMRecord(header);
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
			sw.addAlignment(rec);
			}
		}
	
	private void createEmptyFile(
			SAMFileHeader header,
			String groupName
			) throws IOException
		{
		File fileout=new File(this.OUT_FILE_PATTERN.replaceAll(REPLACE_CHROM, groupName));
		info("creating mock BAM file "+fileout);
		File parent=fileout.getParentFile();
		if(parent!=null) parent.mkdirs();
	
		SAMFileWriter sw=this.samFileWriterFactory.makeBAMWriter(header,false, fileout);
		
		
		if(this.ADD_MOCK_RECORD)
			{
			addMockPair(sw,header);
			}
		sw.close();
		}
	
	private static class ManyToMany
		{
		private java.util.Map<String,Set<String>> group2chroms=new java.util.HashMap<String, java.util.Set<String>>();
		private java.util.Map<String,String> chrom2group=new java.util.HashMap<String, String>();
		public void set(String group,String chrom)
			{
			if(containsChrom(chrom)) throw new IllegalArgumentException("chrom "+chrom+" already defined for group "+ chrom2group.get(chrom));
			java.util.Set<String> set=group2chroms.get(group);
			if(set==null)
				{
				set=new java.util.LinkedHashSet<String>();
				group2chroms.put(group,set);
				}
			set.add(chrom);
			chrom2group.put(chrom,group);
			}
		public boolean containsGroup(String s)
			{
			return group2chroms.containsKey(s);
			}
		public boolean containsChrom(String s)
			{
			return chrom2group.containsKey(s);
			}
			
		}
	
	private void scan(SamReader samFileReader,File chromGroupFile) throws Exception
		{
		SAMSequenceDictionary samSequenceDictionary=samFileReader.getFileHeader().getSequenceDictionary();
		if(samSequenceDictionary==null || samSequenceDictionary.isEmpty())
			{
			throw new PicardException("input is missing a sequence dictionary");
			}
		ManyToMany many2many=new ManyToMany();

		many2many.set(this.UNDERTERMINED_NAME, this.UNDERTERMINED_NAME);

		
		
		if(chromGroupFile!=null)
			{
			Set<String> all_chromosomes=new HashSet<String>();

			for(SAMSequenceRecord seq:samSequenceDictionary.getSequences())
				{
				all_chromosomes.add(seq.getSequenceName());
				}
			
			BufferedReader r=IOUtils.openFileForBufferedReading(chromGroupFile);
			String line;
			while((line=r.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[] =line.split("[ \t,]+");
				String groupName=tokens[0].trim();
				if(groupName.isEmpty()) throw new IOException("Empty group name in "+line);
				if(many2many.containsGroup(groupName))  throw new IOException("Group defined twice "+groupName);
				for(int i=1;i< tokens.length;i++)
					{
					String chromName=tokens[i].trim();
					if(!all_chromosomes.contains(chromName))
						{
						 throw new IOException("chrom "+chromName+" undefined in ref dict");
						}
					if(chromName.isEmpty()) continue;
					many2many.set(groupName,chromName);
					}
				}
			r.close();
			}
		
		for(SAMSequenceRecord seq:samSequenceDictionary.getSequences())
			{
			String chromName=seq.getSequenceName();
			if(many2many.containsChrom(chromName)) continue;
			if(many2many.containsGroup(chromName)) 
				{
				throw new IOException("cannot create chrom group "+chromName+" because it is already defined.");
				}
			many2many.set(chromName,chromName);
			}
		
		
		Map<String,SAMFileWriter> seen=new HashMap<String,SAMFileWriter>(many2many.group2chroms.size());
		
		final SAMFileHeader header=samFileReader.getFileHeader();
		if(createIndex)
			{
			header.setSortOrder(SortOrder.coordinate);
			this.samFileWriterFactory.setCreateIndex(true);
			}
		/*
		 problem of parsing with GATK 2.6 : ignore this for the moment.
		SAMProgramRecord sp=new SAMProgramRecord(getClass().getSimpleName());
		sp.setProgramName(getClass().getSimpleName());
		sp.setProgramVersion(String.valueOf(getProgramVersion()));
		sp.setPreviousProgramGroupId(getClass().getSimpleName());
		sp.setCommandLine(getCommandLine().replaceAll("[ \\s]+"," "));
		header.addProgramRecord(sp);
*/
		
		
        
        
       
       SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(samFileReader.getFileHeader()==null?null:samFileReader.getFileHeader().getSequenceDictionary());
        
		for(Iterator<SAMRecord> iter=samFileReader.iterator();
				iter.hasNext(); )
			{
			SAMRecord record=iter.next();
			progress.watch(record);
			String recordChromName=null;
			if( record.getReadUnmappedFlag() )
				{
				if(!record.getReadPairedFlag() || record.getMateUnmappedFlag())
					{
					recordChromName=this.UNDERTERMINED_NAME;
					}
				else
					{
					recordChromName=record.getMateReferenceName();
					}
				}
			else
				{
				recordChromName=record.getReferenceName();
				
				}
			String groupName=many2many.chrom2group.get(recordChromName);
			if(groupName==null)
				{
				samFileReader.close();
				throw new IOException("Undefined group/chrom for "+recordChromName+" (not in ref dictionary "+many2many.chrom2group.keySet()+").");
				}
			
			SAMFileWriter writer=seen.get(groupName);
			
			if(writer==null)
				{
				File fileout=new File(this.OUT_FILE_PATTERN.replaceAll(REPLACE_CHROM, groupName));
				info("opening "+fileout);
				File parent=fileout.getParentFile();
				if(parent!=null) parent.mkdirs();

				writer=this.samFileWriterFactory.makeBAMWriter(
						header,
						false,
						fileout
						);
				seen.put(groupName, writer);
				}
			writer.addAlignment(record);
			}
		progress.finish();
		
		for(String k:seen.keySet())
			{
			info("closing group "+k);
			seen.get(k).close();
			}
		samFileReader.close();
		
		if(this.GENERATE_EMPTY_BAM)
			{
			
			for(String groupName:many2many.group2chroms.keySet())
				{
				if(seen.containsKey(groupName)) continue;
				createEmptyFile(header,groupName);
				}
			}
		}
	
	@Override
	public String getProgramDescription() {
		return "Split a BAM by chromosome group.";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/SplitBam";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-p (file) output file pattern. MUST contain "+REPLACE_CHROM+" and end with .bam");
		out.println("-g (file) Chromosome group file. Optional");
		out.println("-T (dir) add tmp directory. Optional");
		out.println("-u (name) unmapped chromosome name. Optional. Default "+UNDERTERMINED_NAME);
		out.println("-m add mock record if no samRecord saved in bam");
		out.println("-E generate empty bam if no samRecord found for a given group.");
		out.println("-S sort/create index");
		out.println("-R (int)  max records in RAM "+getMessageBundle("max.records.in.ram"));
		super.printOptions(out);
		}
	private boolean createIndex=false;

	@Override
	public int doWork(String[] args)
		{
		int maxRecordsInRam=1000000;
		File chromGroupFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"u:g:p:mET:SR:"))!=-1)
			{
			switch(c)
				{	
				case 'R': maxRecordsInRam= Integer.parseInt(opt.getOptArg());break;
				case 'u': UNDERTERMINED_NAME= opt.getOptArg();break;
				case 'p': OUT_FILE_PATTERN= opt.getOptArg();break;
				case 'g': chromGroupFile=new File(opt.getOptArg());break;
				case 'T': addTmpDirectory(new File(opt.getOptArg()));break;
				case 'm': ADD_MOCK_RECORD=true;this.GENERATE_EMPTY_BAM=true;break;
				case 'E': GENERATE_EMPTY_BAM=true;break;
				case 'S': createIndex=true;break;
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
		samFileWriterFactory.setCreateIndex(this.createIndex); 
		SamReader sfr=null;
		try
			{
			
			
			if(!OUT_FILE_PATTERN.contains(REPLACE_CHROM))
				{
				error("output file pattern undefined or doesn't contain "+REPLACE_CHROM);
				return -1;
				}
			if(!OUT_FILE_PATTERN.endsWith(".bam"))
				{
				error("output file must end with '.bam'");
				return -1;
				}
			
			samFileWriterFactory.setTempDirectory(  this.getTmpDirectories().get(0) );
			samFileWriterFactory.setMaxRecordsInRam(maxRecordsInRam);
			
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
		new SplitBam2().instanceMainWithExit(args);
		}
	
	}
