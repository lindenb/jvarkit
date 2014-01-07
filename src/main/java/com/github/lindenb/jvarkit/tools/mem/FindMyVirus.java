package com.github.lindenb.jvarkit.tools.mem;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import net.sf.picard.PicardException;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;
import net.sf.samtools.SAMFileWriter;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFlag;
import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;
import com.github.lindenb.jvarkit.util.picard.XPAlign;
import com.github.lindenb.jvarkit.util.picard.XPalignFactory;

public class FindMyVirus extends AbstractCommandLineProgram
	{
	/* how to split the bam, witch categories */
	private enum CAT{
		both_ref()
			{
			@Override
			public String getDescription() {
				return "Both reads map the HOST reference";
				}
			} ,
		both_virus()
			{
			@Override
			public String getDescription() {
				return "Both reads map the VIRUS reference";
				}
			},
		ref_orphan()
			{
			@Override
			public String getDescription() {
				return "In a pair: the read is mapped on the HOST Reference, the mate is unmapped. No spliced read mapped on VIRUS.";
				}
			},
		virus_orphan()
			{
			@Override
			public String getDescription() {
				return "In a pair: the read is mapped on the VIRUS Reference, the mate is unmapped. No spliced read mapped on HOST.";
				}
			},
		ref_and_virus()
			{
			@Override
			public String getDescription() {
				return "In a pair: a read is mapped on the VIRUS Reference, the other is mapped on HOST. No spliced read mapped on the other genome.";
				}
			},
		ref_and_virus_spliced()
			{
			@Override
			public String getDescription() {
				return "In a pair: a read is mapped on the VIRUS Reference, the other is mapped on HOST. A spliced read is mapped on the other genome.";
				}
			},
		unpaired()
			{
			@Override
			public String getDescription() {
				return "read is not paired.";
				}
			},
		undetermined()
			{
			@Override
			public String getDescription() {
				return "Unknown category/case";
				}
			},
		unmapped()
			{
			@Override
			public String getDescription() {
				return "Read is unmapped";
				}
			},
		duplicate()
			{
			@Override
			public String getDescription() {
				return "Read is duplicate";
				}
			},
		secondary
			{
			@Override
			public String getDescription() {
				return "Secondary alignment";
				}
			},
		failsqual{
			@Override
			public String getDescription() {
				return "Fails mapper quality-control";
				}
			};
		
		public String getDescription()
			{
			return name();
			}
		};
	@Override
	public String getProgramDescription() {
		return "Find my Virus. Created for Adrien Inserm. Proportion of reads mapped on HOST/VIRUS. ";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -o (file) bam out prefix");
		out.println(" -V (name) virus chrom name (may be called multiple times)");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		Set<String> virusNames=new HashSet<String>();
		File bamOut=null;
		SamWriterFactory sfwf= SamWriterFactory.newInstance();
		sfwf.setBinary(true);
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:V:"))!=-1)
			{
			switch(c)
				{
				case 'V': virusNames.add(opt.getOptArg());break;
				case 'o': bamOut=new File(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(virusNames.isEmpty())
			{
			error("no virus name");
			return -1;
			}
		if(bamOut==null || (bamOut.exists() && bamOut.isDirectory()))
			{
			error("Bam file undefined");
			return -1;
			}
		
		
		
		SAMFileReader sfr=null;
		SAMFileWriter sfwArray[]=new SAMFileWriter[CAT.values().length];
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				sfr=new SAMFileReader(new File(args[opt.getOptInd()]));
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=sfr.getFileHeader();
			for(CAT category:CAT.values())
				{
				info("Opening "+category);
				SAMFileHeader header2=header.clone();
				header2.addComment("Category:"+category.name());
				header2.addComment("Description:"+category.getDescription());
				SAMProgramRecord rec=header2.createProgramRecord();
				rec.setCommandLine(this.getProgramCommandLine());
				rec.setProgramName(getProgramName());
				rec.setProgramVersion(getVersion());
				rec.setAttribute("CAT", category.name());
				File outputFile=new File(bamOut.getParentFile(),bamOut.getName()+"."+category.name()+".bam");
				info("Opening "+outputFile);
				File countFile=new File(bamOut.getParentFile(),bamOut.getName()+"."+category.name()+".count.txt");
				SAMFileWriter sfw=sfwf.make(header2, outputFile);
				sfw=new SAMFileWriterCount(sfw, countFile,category);
				sfwArray[category.ordinal()]=sfw;
				}
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			XPalignFactory xpAlignFactory=new XPalignFactory(header);
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				progress.watch(rec);
				CAT category=null;
				List<XPAlign> xpList=Collections.emptyList();
				if(category==null && !rec.getReadPairedFlag())
					{
					category=CAT.unpaired;
					}
				
				if(category==null && rec.isSecondaryOrSupplementary())
					{
					category=CAT.secondary;
					}
				
				if(category==null &&  rec.getReadFailsVendorQualityCheckFlag())
					{
					category=CAT.failsqual;
					}
				if(category==null &&  rec.getDuplicateReadFlag())
					{	
					category=CAT.duplicate;
					}

				
				
				if(category==null &&  rec.getReadUnmappedFlag())
					{
					category=CAT.unmapped;
					}
				
				if(category==null)
					{
					xpList=xpAlignFactory.getXPAligns(rec);
					}
			
				boolean xp_containsVirus=false;
				boolean xp_containsChrom=false;
				for(XPAlign xpa:xpList)
					{
					if(virusNames.contains(xpa.getChrom()))
						{
						xp_containsVirus=true;
						}
					else
						{
						xp_containsChrom=true;
						}
					}
				
				
				/* both reads mapped on ref */
				if(category==null &&
					!rec.getReadUnmappedFlag() &&
					!rec.getMateUnmappedFlag() &&
					!virusNames.contains(rec.getReferenceName()) &&
					!virusNames.contains(rec.getMateReferenceName())
					)
					{
					if(!xp_containsVirus)
						{
						category=CAT.both_ref;
						}
					else
						{
						category=CAT.ref_and_virus_spliced;
						}
					}
				
				/*  pair(unmapped,mapped on reference) */
				if(category==null &&
						(
						 (!rec.getReadUnmappedFlag() && rec.getMateUnmappedFlag() && !virusNames.contains(rec.getReferenceName())) ||
						 (rec.getReadUnmappedFlag() && !rec.getMateUnmappedFlag() && !virusNames.contains(rec.getMateReferenceName())) 
						)) 
						{
						if(!xp_containsVirus)
							{
							category=CAT.ref_orphan;
							}
						else
							{
							category=CAT.ref_and_virus_spliced;
							}
						}
				
				
				/* both reads mapped on virus */
				if(category==null &&
					!rec.getReadUnmappedFlag() &&
					!rec.getMateUnmappedFlag() &&
					virusNames.contains(rec.getReferenceName()) &&
					virusNames.contains(rec.getMateReferenceName())
					)
					{
					if(!xp_containsChrom)
						{
						category=CAT.both_virus;
						}
					else
						{
						category=CAT.ref_and_virus_spliced;
						}
					}

				
				if(category==null &&
						(
						 (!rec.getReadUnmappedFlag() && rec.getMateUnmappedFlag() && virusNames.contains(rec.getReferenceName())) ||
						 (rec.getReadUnmappedFlag() && !rec.getMateUnmappedFlag() && virusNames.contains(rec.getMateReferenceName())) 
						)) 
						{
						if(!xp_containsChrom)
							{
							category=CAT.virus_orphan;
							}
						else
							{
							category=CAT.ref_and_virus_spliced;
							}
						}
				
				
				
				if(category==null &&
					!rec.getReadUnmappedFlag() &&
					!rec.getMateUnmappedFlag() &&
					(
					(virusNames.contains(rec.getReferenceName()) && !virusNames.contains(rec.getMateReferenceName())) ||
					(!virusNames.contains(rec.getReferenceName()) && virusNames.contains(rec.getMateReferenceName()))
					)
					)
					{
					category=CAT.ref_and_virus;
					}
				
				
				
				/*dispatch */
				if(category==null)
					{
					warning("Not handled: "+rec);
					category=CAT.undetermined;
					}
				sfwArray[category.ordinal()].addAlignment(rec);
				}
			for(SAMFileWriter sfw:sfwArray)
				{
				info("Closing "+sfw);
				sfw.close();
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
			info("Closing");
			CloserUtil.close(sfr);
			CloserUtil.close(sfwArray);
			}
		}
	
	public static void main(String[] args)
		{
		new FindMyVirus().instanceMainWithExit(args);
		}
	
	/** wrapper of SAMFileWriter to get some statistics about a BAM. write data on close() */
	private class SAMFileWriterCount
		implements SAMFileWriter
		{
		private CAT category;
		private SAMFileWriter delegate;
		private File countFile;
		private Counter<String> chrom=new Counter<String>();
		private Counter<Integer> flags=new Counter<Integer>();
		SAMFileWriterCount(SAMFileWriter delegate,File countFile,CAT category)
			{
			this.category=category;
			this.countFile=countFile;
			this.delegate=delegate;
			for(SAMSequenceRecord rec:delegate.getFileHeader().getSequenceDictionary().getSequences())
				{
				chrom.initializeIfNotExists(rec.getSequenceName());
				}
			chrom.initializeIfNotExists(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
			}
		@Override
		public void addAlignment(SAMRecord alignment) {
			this.delegate.addAlignment(alignment);
			this.flags.incr(alignment.getFlags());
			if(!alignment.getReadUnmappedFlag())
				{
				chrom.incr(alignment.getReferenceName());
				}
			else
				{
				chrom.incr(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				}
			}
		@Override
		public SAMFileHeader getFileHeader() {
			return this.delegate.getFileHeader();
			}
		@Override
		public void close()
			{
			info("Closing SAMFileWriterCount ");
			
			PrintWriter fw=null;
			try
				{
				info("Writing "+countFile);
				fw=new PrintWriter(countFile);
				fw.println(this.category.name());
				fw.println(this.category.getDescription());
				fw.println("#CHROMOSOME\tCOUNT");
				for(String c:this.chrom.keySetDecreasing())
					{
					fw.println(c+"\t"+this.chrom.count(c));
					}
				fw.println("Total\t"+this.chrom.getTotal());
				
				fw.println("#FLAG\tCOUNT\texplain");
				for(Integer c:this.flags.keySetDecreasing())
					{
					fw.print(c+"\t"+this.flags.count(c)+"\t");
					for(SamFlag flg:SamFlag.values())
						{
						if(flg.isSet(c))
							{
							fw.write(flg.name()+" ");
							}
						}
					fw.println();
					}
				fw.flush();
				}
			catch(Exception err)
				{
				error(err);
				throw new PicardException("Boum:"+countFile,err);
				}
			finally
				{
				CloserUtil.close(fw);
				}
			this.delegate.close();
			}
		@Override
		public String toString() {
			return "SAMFileWriterCount "+countFile;
			}
		}

	}
