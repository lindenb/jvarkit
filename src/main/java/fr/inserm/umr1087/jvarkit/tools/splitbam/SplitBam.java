package fr.inserm.umr1087.jvarkit.tools.splitbam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.DefaultSAMRecordFactory;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecordFactory;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class SplitBam
	{
	private static final Logger LOG=Logger.getLogger("split");
	private final String REPLACE_CHROM="__CHROM__";
	private String outFilePattern="";
	private String underterminedName="Unmapped";
	private boolean generate_empty_bams=false;
	private SAMSequenceDictionary  samSequenceDictionary;
	private boolean if_bam_empty_add_mock_sam_record=false;
	private long id_generator=System.currentTimeMillis();
	private boolean input_is_sorted=false;
	private boolean create_index=false;
	private File tmpDir=null;
	private File chromGroup=null;
	
	
	
	private SplitBam()
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
			rec.setReadBases(bases.getBytes());
			rec.setMappingQuality(0);
			rec.setBaseQualityString(bases.replace('N', 'I'));
			rec.setReadUnmappedFlag(true);
			rec.setMateUnmappedFlag(true);
			rec.setReadPairedFlag(true);
			rec.setReadName("MOCKREAD"+(id_generator)+":6:190:289:82");
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
		File fileout=new File(this.outFilePattern.replaceAll(REPLACE_CHROM, groupName));
		LOG.info("creating mock BAM file "+fileout);
		File parent=fileout.getParentFile();
		if(parent!=null) parent.mkdirs();

		SAMFileWriter sw=sf.makeBAMWriter(header, true, fileout);
		if(if_bam_empty_add_mock_sam_record)
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

		many2many.set(this.underterminedName, this.underterminedName);

		
		
		if(this.chromGroup!=null)
			{
			Set<String> all_chromosomes=new HashSet<String>();

			for(SAMSequenceRecord seq:this.samSequenceDictionary.getSequences())
				{
				all_chromosomes.add(seq.getSequenceName());
				}
			
			BufferedReader r=new BufferedReader(new FileReader(this.chromGroup));
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
		SAMFileReader samFileReader=new SAMFileReader(in);
		samFileReader.setValidationStringency(ValidationStringency.SILENT);
		SAMFileHeader header=samFileReader.getFileHeader();
		header.setSortOrder(SortOrder.coordinate);
		
        SAMFileWriterFactory sf=new SAMFileWriterFactory();
       if(this.tmpDir!=null)
    	   {
    	   sf.setTempDirectory(this.tmpDir);
    	   }
        sf.setCreateIndex(this.create_index);
        
        long nrecords=0L;
       
       
        
		for(Iterator<SAMRecord> iter=samFileReader.iterator();
				iter.hasNext(); )
			{
			SAMRecord record=iter.next();
			++nrecords;
			if(nrecords%1E6==0)
				{
				LOG.info("nRecord:"+nrecords);
				}
			String recordChromName=null;
			if( record.getReadUnmappedFlag() )
				{
				if(record.getMateUnmappedFlag())
					{
					recordChromName=this.underterminedName;
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
				throw new IOException("Undefined group/chrom for "+recordChromName+" (not in ref dictionary "+many2many.chrom2group.keySet()+").");
				}
			
			SAMFileWriter writer=seen.get(groupName);
			
			if(writer==null)
				{
				File fileout=new File(this.outFilePattern.replaceAll(REPLACE_CHROM, groupName));
				LOG.info("opening "+fileout);
				File parent=fileout.getParentFile();
				if(parent!=null) parent.mkdirs();
				writer=sf.makeBAMWriter(header,this.input_is_sorted,fileout);
				seen.put(groupName, writer);
				nrecords=0L;
				}
			
			writer.addAlignment(record);
			}
		
		for(String k:seen.keySet())
			{
			LOG.info("closing group "+k);
			seen.get(k).close();
			}
		samFileReader.close();
		
		if(generate_empty_bams)
			{
			
			for(String groupName:many2many.group2chroms.keySet())
				{
				if(seen.containsKey(groupName)) continue;
				createEmptyFile(sf,header,groupName);
				}
			}
		
		}
	
	
	private void run(String[] args)
		throws Exception
		{
		File referenceFile=null;
		
		int optind=0;
		while(optind< args.length)
			{
			if(args[optind].equals("-h") ||
			   args[optind].equals("-help") ||
			   args[optind].equals("--help"))
				{
				System.err.println("Pierre Lindenbaum PhD. 2013");
				System.err.println("Options:");
				System.err.println(" -h help; This screen.");
				System.err.println(" -R (reference file) REQUIRED.");
				System.err.println(" -u (unmapped chromosome name): default:"+this.underterminedName);
				System.err.println(" -e | --empty : generate EMPTY bams for chromosome having no read mapped");
				System.err.println(" -m | --mock : if option '-e', add a mock pair of sam records to the bam");
				System.err.println(" -p (output file/bam pattern) REQUIRED. MUST contain "+REPLACE_CHROM+" and end with .bam");
				System.err.println(" -s assume input is sorted.");
				System.err.println(" -x | --index  create index.");
				System.err.println(" -t | --tmp  (dir) tmp file directory");
				System.err.println(" -G (file) chrom-group file\n" +
							       "     Merge some chromosome in the following groups. Format:\n" +
							       "     (group-name1)\\tchrom1\\tchrom2\\tchrom3...\\n\n"+
							       "     (group-name2)\\tchrom11\\tchrom12\\tchrom22...\\n\n"+
							       "     The missing chromosomes are defined in their own group"
									);
				return;
				}
			else if(args[optind].equals("-e")|| args[optind].equals("--empty"))
				{
				this.generate_empty_bams=true;
				}
			else if(args[optind].equals("-m") || args[optind].equals("--mock"))
				{
				this.generate_empty_bams=true;
				this.if_bam_empty_add_mock_sam_record=true;
				}
			else if((args[optind].equals("-t") || args[optind].equals("--tmp")) && optind+1< args.length)
				{
				this.tmpDir=new File(args[++optind]);
				}
			else if(args[optind].equals("-R") && optind+1< args.length)
				{
				referenceFile=new File(args[++optind]);
				}
			else if(args[optind].equals("-u") && optind+1< args.length)
				{
				underterminedName= args[++optind];
				}
			else if(args[optind].equals("-p") && optind+1< args.length)
				{
				outFilePattern= args[++optind];
				}
			else if(args[optind].equals("-G") && optind+1< args.length)
				{
				chromGroup= new File(args[++optind]);
				}
			else if(args[optind].equals("-s"))
				{
				this.input_is_sorted=true;
				}
			else if(args[optind].equals("-x") || args[optind].equals("--index"))
				{
				this.create_index=true;
				}
			else if(args[optind].equals("--"))
				{
				optind++;
				break;
				}
			else if(args[optind].startsWith("-"))
				{
				System.err.println("Unknown option "+args[optind]);
				return;
				}
			else 
				{
				break;
				}
			++optind;
			}
		if(!outFilePattern.contains(REPLACE_CHROM))
			{
			System.err.println("output file pattern undefined or doesn't contain "+REPLACE_CHROM);
			System.exit(-1);
			}
		
		if(referenceFile==null)
			{
			System.err.println("Reference file undefined");
			System.exit(-1);
			}
		IndexedFastaSequenceFile indexedFastaSequenceFile=new IndexedFastaSequenceFile(referenceFile);
		this.samSequenceDictionary=indexedFastaSequenceFile.getSequenceDictionary();
		if(this.samSequenceDictionary==null)
			{
			System.err.println("Reference file dictionary missing. use picard to create it.");
			System.exit(-1);
			}
		
		if(optind==args.length)
			{
			scan(System.in);
			}
		else if(optind+1==args.length)
			{
			FileInputStream fin=new FileInputStream(args[optind]);
			scan(fin);
			fin.close();
			}
		else 
			{
			System.err.println("illegal number of arguments.");
			System.exit(-1);
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		new SplitBam().run(args);
		}
	
	}
