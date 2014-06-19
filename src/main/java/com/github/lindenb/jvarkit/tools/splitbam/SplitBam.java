package com.github.lindenb.jvarkit.tools.splitbam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryFactory;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;

import com.github.lindenb.jvarkit.util.picard.PicardException;
import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.StandardOptionDefinitions;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.BAMIndex;
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
import htsjdk.samtools.util.SequenceUtil;

/***
 * 
 * SplitBam
 *
 */
@Deprecated
public class SplitBam extends AbstractCommandLineProgram
	{
	
	private static final Log LOG=Log.getInstance(SplitBam.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+"Split a BAM by chromosome group.";

    @Option(shortName= StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Indexex reference",optional=false)
	public File  REF=null;
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="BAM file to process. Default stdin. ",optional=true)
	public File IN=null;
    
    @Option(shortName= "EMPTY_BAM", doc="generate EMPTY bams for chromosome having no read mapped. ",optional=true)
	public boolean GENERATE_EMPTY_BAM=false;

    @Option(shortName= "GP", doc="Chromosome group file. ",optional=true)
	public File CHROM_GROUP=null;
    
    
    @Option(shortName= "MOCK", doc="add a mock pair of sam records to the bam. ",optional=true)
	public boolean ADD_MOCK_RECORD=false;

    
	private final static String REPLACE_CHROM="__CHROM__";
	
	
	@Option(shortName= "OFP", doc="MUST contain "+REPLACE_CHROM+" and end with .bam. ",optional=false)
	public String OUT_FILE_PATTERN="";
	
	@Option(shortName= "UN", doc="Unmapped chromosome name. ",optional=true)
	public String UNDERTERMINED_NAME="Unmapped";
	
	@Option(shortName= "IS", doc="input is sorted. ",optional=true)
	public boolean INPUT_IS_SORTED=false;
	
	private SAMSequenceDictionary  samSequenceDictionary;
	private long id_generator=System.currentTimeMillis();
	
	private Set<File> deleteOnError=new HashSet<File>();

	
	public SplitBam()
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
			String readName="MOCKREAD"+(id_generator)+":6:190:289:82";
			rec.setReadName(readName);
			LOG.info("generating mock read: "+readName);
			rec.setAttribute("MK",1);
			if(G!=null && !G.isEmpty())
				{
				rec.setAttribute("RG", G.get(0).getId());
				}
			sw.addAlignment(rec);
			}
		}
	
	private void createEmptyFile(
			SAMFileWriterFactory sf,
			SAMFileHeader header,
			String groupName
			) throws IOException
		{
		File fileout=new File(this.OUT_FILE_PATTERN.replaceAll(REPLACE_CHROM, groupName));
		LOG.info("creating mock BAM file "+fileout);
		File parent=fileout.getParentFile();
		if(parent!=null) parent.mkdirs();

		SAMFileWriter sw=sf.makeBAMWriter(header, true, fileout);
		
		deleteOnError.add(fileout);
		deleteOnError.add(indexFor(fileout));
		
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
	
	private void scan(InputStream in) throws Exception
		{
		ManyToMany many2many=new ManyToMany();

		many2many.set(this.UNDERTERMINED_NAME, this.UNDERTERMINED_NAME);

		
		
		if(this.CHROM_GROUP!=null)
			{
			Set<String> all_chromosomes=new HashSet<String>();

			for(SAMSequenceRecord seq:this.samSequenceDictionary.getSequences())
				{
				all_chromosomes.add(seq.getSequenceName());
				}
			
			BufferedReader r=IOUtils.openFileForBufferedReading(this.CHROM_GROUP);
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
		
		for(SAMSequenceRecord seq:this.samSequenceDictionary.getSequences())
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
		SamReader samFileReader=SamFileReaderFactory.mewInstance().open(in);
		final SAMFileHeader header=samFileReader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		
		if(!SequenceUtil.areSequenceDictionariesEqual(
				header.getSequenceDictionary(),
				this.samSequenceDictionary)
				)
			{
			samFileReader.close();
			throw new PicardException("Not the same sequence dictionary BAM vs "+REF);
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
		
		
        SAMFileWriterFactory sf=new SAMFileWriterFactory();
       if(!super.TMP_DIR.isEmpty())
    	   {
    	   sf.setTempDirectory(super.TMP_DIR.get(0));
    	   }
        if(super.CREATE_INDEX!=null)
        	{
        	sf.setCreateIndex(super.CREATE_INDEX);
        	}
        
       
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
				LOG.info("opening "+fileout);
				File parent=fileout.getParentFile();
				if(parent!=null) parent.mkdirs();
				writer=sf.makeBAMWriter(
						header,
						this.INPUT_IS_SORTED,
						fileout,
						super.COMPRESSION_LEVEL
						);
				deleteOnError.add(fileout);
				deleteOnError.add(indexFor(fileout));

				
				seen.put(groupName, writer);
				
				
				
				}
			writer.addAlignment(record);
			}
		
		for(String k:seen.keySet())
			{
			LOG.info("closing group "+k);
			seen.get(k).close();
			}
		samFileReader.close();
		
		if(this.GENERATE_EMPTY_BAM)
			{
			
			for(String groupName:many2many.group2chroms.keySet())
				{
				if(seen.containsKey(groupName)) continue;
				createEmptyFile(sf,header,groupName);
				}
			}
		deleteOnError.clear();
		}
	@Override
	protected int doWork()
		{
		try
			{
			
			
			
			if(this.ADD_MOCK_RECORD)
				{
				this.GENERATE_EMPTY_BAM=true;
				}
			
			if(!OUT_FILE_PATTERN.contains(REPLACE_CHROM))
				{
				LOG.error("output file pattern undefined or doesn't contain "+REPLACE_CHROM);
				return -1;
				}
			
			if(REF==null)
				{
				LOG.error("Reference file undefined");
				System.exit(-1);
				}
			this.samSequenceDictionary=new SAMSequenceDictionaryFactory().load(REF);
			if(this.samSequenceDictionary==null)
				{
				LOG.error("Reference file dictionary missing. use picard to create it.");
				return -1;
				}
			
			if(this.IN==null)
				{
				LOG.info("reading stdin");
				scan(System.in);
				}
			else 
				{
				LOG.info("reading "+IN);
				FileInputStream fin=new FileInputStream(IN);
				scan(fin);
				fin.close();
				}
			return 0;
			}
		catch(Exception err)
			{
			err.printStackTrace();
			for(File f:deleteOnError)
				{
				if(f!=null && f.exists() && f.isFile())
					{
					LOG.info("Cleanup: delete "+f);
					f.delete();
					}
				}
			super.testRemoteGit();
			return -1;
			}
		finally
			{
			
			}
		}
	
	private static File indexFor(File bamFile)
		{
		final String bamExtension = ".bam";
		String fileName=bamFile.getName();
		final String bai = fileName.substring(0, fileName.length() - bamExtension.length()) + BAMIndex.BAMIndexSuffix;
        return new File(bamFile.getParent(), bai);
		}
	
	@Override
	public String getProgramVersion() {
		return "1.0";
		}
	
	public static void main(String[] args) throws Exception
		{
		new SplitBam().instanceMainWithExit(args);
		}
	
	}
