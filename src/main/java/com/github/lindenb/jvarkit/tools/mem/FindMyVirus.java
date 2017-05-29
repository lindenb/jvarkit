package com.github.lindenb.jvarkit.tools.mem;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLoggerInterface;
import htsjdk.samtools.SAMFileWriter;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;

/**

## Example

```bash
$  java -jar dist/findmyvirus.jar -V virus_chr -o category in.bam
```


*/
@Program(name="findmyvirus",
	description="Find my Virus. Created for Adrien Inserm. Proportion of reads mapped on HOST/VIRUS. ")
public class FindMyVirus extends Launcher
	{
	private static final Logger LOG= Logger.build(FindMyVirus.class).make();
	
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

	@Parameter(names={"-o","--out"},description="output",required=true)
	private File outputFile=null;

	@Parameter(names={"-V"},description=" virus chrom name")
	private	Set<String> virusNames=new HashSet<String>();

	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	@Override
	public int doWork(List<String> args) {
		
		if(virusNames.isEmpty())
			{
			LOG.error("no virus name");
			return -1;
			}
		
		
		
		
		SamReader sfr=null;
		SAMFileWriter sfwArray[]=new SAMFileWriter[CAT.values().length];
		try
			{
			sfr = openSamReader(oneFileOrNull(args));
			SAMFileHeader header=sfr.getFileHeader();
			for(CAT category:CAT.values())
				{
				LOG.info("Opening "+category);
				SAMFileHeader header2=header.clone();
				header2.addComment("Category:"+category.name());
				header2.addComment("Description:"+category.getDescription());
				SAMProgramRecord rec=header2.createProgramRecord();
				rec.setCommandLine(this.getProgramCommandLine());
				rec.setProgramName(getProgramName());
				rec.setProgramVersion(getVersion());
				rec.setAttribute("CAT", category.name());
				File outputFile=new File(this.outputFile.getParentFile(),this.outputFile.getName()+"."+category.name()+".bam");
				LOG.info("Opening "+outputFile);
				File countFile=new File(outputFile.getParentFile(),outputFile.getName()+"."+category.name()+".count.txt");
				SAMFileWriter sfw= writingBamArgs.openSAMFileWriter(outputFile, header2,  true);
				sfw=new SAMFileWriterCount(sfw, countFile,category);
				sfwArray[category.ordinal()]=sfw;
				}
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			OtherCanonicalAlignFactory xpAlignFactory=new OtherCanonicalAlignFactory(header);
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				progress.watch(rec);
				CAT category=null;
				List<OtherCanonicalAlign> xpList=Collections.emptyList();
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
				for(OtherCanonicalAlign xpa:xpList)
					{
					if(virusNames.contains(xpa.getReferenceName()))
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
					LOG.warning("Not handled: "+rec);
					category=CAT.undetermined;
					}
				sfwArray[category.ordinal()].addAlignment(rec);
				}
			for(SAMFileWriter sfw:sfwArray)
				{
				LOG.info("Closing "+sfw);
				sfw.close();
				}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			LOG.info("Closing");
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
		@SuppressWarnings("unused")
		private ProgressLoggerInterface progressLogger;
		
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
		public void setProgressLogger(ProgressLoggerInterface progressLogger) {
			this.progressLogger=progressLogger;
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
			LOG.info("Closing SAMFileWriterCount ");
			
			PrintWriter fw=null;
			try
				{
				LOG.info("Writing "+countFile);
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
					for(SAMFlag flg:SAMFlag.values())
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
				LOG.error(err);
				throw new RuntimeException("Boum:"+countFile,err);
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
