package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import net.sf.samtools.SAMFileHeader.SortOrder;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloserUtil;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class ConvertBamChromosomes
	extends AbstractCommandLineProgram
	{
	private boolean use_original_chrom_name_if_no_mapping=false;
	private Map<String,String> customMapping=new HashMap<String,String>();
	private Set<String> unmappedChromosomes=new HashSet<String>();
	private boolean ignore_if_no_mapping=false;
	private ConvertBamChromosomes()
		{
		
		}
	
	private String convertName(String chrom)throws IOException
		{
		if(chrom==null) throw new NullPointerException();
		String newname=customMapping.get(chrom);
		if(newname==null)
			{
			if(!unmappedChromosomes.contains(chrom))
				{
				warning("unmapped chromosome "+chrom);
				unmappedChromosomes.add(chrom);
				}
			if(ignore_if_no_mapping) return null;
			
			if(use_original_chrom_name_if_no_mapping)
				{	
				return chrom;
				}
			throw new IOException("No mapping found to convert name of chromosome \""+chrom+"\"");
			}
		return newname;
		}
	
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BamRenameChromosomes";
		}
	
	@Override
	public String getProgramDescription() {
		return "Convert the names of the chromosomes in a BAM file.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -f (file) load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+");
		out.println(" -i if no mapping found, skip that record.");
		out.println(" -C if no mapping found, use the original name instead of throwing an error. ");
		out.println(" -o (filenameout.bam) default: SAM as stdout");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File bamout=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"f:Cio:"))!=-1)
			{
			switch(c)
				{
				case 'o': bamout=new File(opt.getOptArg());break;
				case 'i': ignore_if_no_mapping=true;break;
				case 'C': use_original_chrom_name_if_no_mapping=true;break;
				case 'f':
					{
					File f=new File(opt.getOptArg());
					BufferedReader in=null;
					try
						{
						info("Loading custom mapping "+f);
						in=IOUtils.openFileForBufferedReading(f);
						String line;
						while((line=in.readLine())!=null)
							{
							if(line.isEmpty() || line.startsWith("#")) continue;
							String tokens[]=line.split("[\t]");
							if(tokens.length!=2
									|| tokens[0].trim().isEmpty()
									|| tokens[1].trim().isEmpty()
									|| tokens[0].equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
									|| tokens[1].equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)
									) throw new IOException("Bad mapping line: \""+line+"\"");
							tokens[0]=tokens[0].trim();
							tokens[1]=tokens[1].trim();
							if(customMapping.containsKey(tokens[0]))
								{
								throw new IOException("Mapping defined twice for: \""+tokens[0]+"\"");
								}
							customMapping.put(tokens[0], tokens[1]);
							}
						}
					catch(Exception err)
						{
						error(err);
						return -1;
						}
					finally
						{
						CloserUtil.close(in);
						}
					break;
					}
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
		SAMFileReader sfr=null;
		SAMFileWriter sfw=null;
		try
			{
			
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File fin=new File(args[opt.getOptInd()]);
				info("Reading from "+fin);
				sfr=new SAMFileReader(fin);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			
			sfr.setValidationStringency(ValidationStringency.LENIENT);
			SAMFileHeader header1=sfr.getFileHeader();
			if(header1==null)
				{
				error("File header missing");
				return -1;
				}
			
			
			SAMFileHeader header2=header1.clone();
			
			//create new sequence dict
			final SAMSequenceDictionary dict1=header1.getSequenceDictionary();
			if(dict1==null)
				{
				error("Sequence dict missing");
				return -1;
				}
			List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>(dict1.size());
			for(int i=0;i< dict1.size();++i)
				{
				SAMSequenceRecord ssr=dict1.getSequence(i);
				String newName=convertName(ssr.getSequenceName());
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
			SAMProgramRecord prog=header2.createProgramRecord();
			prog.setCommandLine(this.getProgramCommandLine());
			prog.setProgramName(getProgramName());
			prog.setProgramVersion(getVersion());
			
			
			
			
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(dict1);

			SAMFileWriterFactory sfwf=new SAMFileWriterFactory();
			boolean presorted=(header1.getSortOrder()!=null && (header1.getSortOrder()==SortOrder.coordinate || header1.getSortOrder()==SortOrder.queryname));
			if(bamout!=null)
				{
				info("saving to "+bamout);
				sfw=sfwf.makeSAMOrBAMWriter(header2, presorted, bamout);
				}
			else
				{
				sfw=sfwf.makeSAMWriter(header2, presorted, System.out);
				}
			
			
			long num_ignored=0L;
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				SAMRecord rec1=iter.next();
				progress.watch(rec1);
				String newName1=null;
				String newName2=null;
				if(!rec1.getReadUnmappedFlag())
					{
					newName1=convertName(rec1.getReferenceName());
					}
				if(rec1.getReadPairedFlag() && !rec1.getMateUnmappedFlag())
					{
					newName2=convertName(rec1.getMateReferenceName());
					}
				rec1.setHeader(header2);

				if(!rec1.getReadUnmappedFlag())
					{
					if(newName1==null)
						{
						++num_ignored;
						continue;
						}
					rec1.setReferenceName(newName1);
					}
				if(rec1.getReadPairedFlag() && !rec1.getMateUnmappedFlag())
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
				warning("Unmapped chromosomes: "+unmappedChromosomes);
				}
			warning("num ignored read:"+num_ignored);
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
