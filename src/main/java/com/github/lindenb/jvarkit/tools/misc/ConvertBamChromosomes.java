package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

@Program(name="bamrenamechr",description="Convert the names of the chromosomes in a BAM file")
public class ConvertBamChromosomes
	extends Launcher
	{
	private static final Logger LOG = Logger.build(ConvertBamChromosomes.class).make();
	
	
	@Parameter(names={"-c","-convert"},description="What should I do when  a converstion is not found")
	private ContigNameConverter.OnNotFound onNotFound=ContigNameConverter.OnNotFound.RAISE_EXCEPTION;
	@Parameter(names={"-f","--mapping","-m"},description="load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+",required=true)
	private File mappingFile=null;
	@Parameter(names={"-o","--out"},description="output vcf")
	private File output= null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();

	private ConvertBamChromosomes()
		{
		
		}
	
	@Override
	public int doWork(List<String> args) {
		SamReader sfr=null;
		SAMFileWriter sfw=null;
		final Set<String> unmappedChromosomes = new HashSet<>();
		try
			{
			final ContigNameConverter customMapping = ContigNameConverter.fromFile(mappingFile);
			customMapping.setOnNotFound(this.onNotFound);

			
			sfr = super.openSamReader(oneFileOrNull(args));
			
			SAMFileHeader header1=sfr.getFileHeader();
			if(header1==null)
				{
				LOG.error("File header missing");
				return -1;
				}
			
			
			SAMFileHeader header2=header1.clone();
			
			//create new sequence dict
			final SAMSequenceDictionary dict1=header1.getSequenceDictionary();
			if(dict1==null)
				{
				LOG.error("Sequence dict missing");
				return -1;
				}
			final List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>(dict1.size());
			for(int i=0;i< dict1.size();++i)
				{
				SAMSequenceRecord ssr=dict1.getSequence(i);
				String newName=customMapping.apply(ssr.getSequenceName());
				if(newName==null)
					{
					//skip unknown chromosomes
					continue;
					}
				ssr=new SAMSequenceRecord(newName, ssr.getSequenceLength());
				ssrs.add(ssr);
				}
			header2.setSequenceDictionary(new SAMSequenceDictionary(ssrs));
			
			SAMSequenceDictionary dict2=new SAMSequenceDictionary(ssrs);
			header2.setSequenceDictionary(dict2);
			
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1);

			sfw = this.writingBamArgs.openSAMFileWriter(output, header2, true);
			
			
			long num_ignored=0L;
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec1=iter.next();
				progress.watch(rec1);
				String newName1=null;
				String newName2=null;
				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					newName1=customMapping.apply(rec1.getReferenceName());
					}
				if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
					{
					newName2=customMapping.apply(rec1.getMateReferenceName());
					}
				rec1.setHeader(header2);

				if(!SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getReferenceName()))
					{
					if(newName1==null)
						{
						++num_ignored;
						continue;
						}
					rec1.setReferenceName(newName1);
					}
				if(rec1.getReadPairedFlag() && !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(rec1.getMateReferenceName()))
					{
					if(newName2==null)
						{
						++num_ignored;
						continue;
						}
					rec1.setMateReferenceName(newName2);
					}
				sfw.addAlignment(rec1);
				}
			if(!unmappedChromosomes.isEmpty())
				{
				LOG.warning("Unmapped chromosomes: "+unmappedChromosomes);
				}
			LOG.warning("num ignored read:"+num_ignored);
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfr);
			CloserUtil.close(sfw);
			}
		}
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new ConvertBamChromosomes().instanceMainWithExit(args);
		}
	}
