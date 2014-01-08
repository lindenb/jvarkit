package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;

import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;

public class Biostar90204 extends AbstractCommandLineProgram
	{
	private String prefix="_splitbam";
	private int suffix_length=2;
	private  Biostar90204() {
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/Biostar90204";
		}
	
	@Override
	public String getProgramDescription() {
		return "Bam version of linux split. See also http://www.biostars.org/p/90204/";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -n (int) put NUMBER of SAM RECORD per output file. REQUIRED.");
		out.println(" -p (prefix) output file prefix. default:"+this.prefix);
		out.println(" -a (int) suffix length. default:."+suffix_length);
		out.println(" -M (fileout) manifest file. Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		SamWriterFactory swfactory=SamWriterFactory.newInstance();
		File manifestFile=null;
		long record_per_file=-1L;
		com.github.lindenb.jvarkit.util.cli.GetOpt getopt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=getopt.getopt(args,getGetOptDefault()+ "n:p:m:a:"))!=-1)
			{
			switch(c)
				{
				case 'a':
					{
					suffix_length=Integer.parseInt(getopt.getOptArg());
					if(suffix_length<0)
						{
						error("Bad value of suffix_length:"+suffix_length);
						return -1;
						}
					break;
					}
				case 'n':
					{
					record_per_file=Long.parseLong(getopt.getOptArg());
					if(record_per_file<1L)
						{
						error("Bad value of record_per_file:"+record_per_file);
						return -1;
						}
					break;
					}
				case 'p':
					{
					prefix=getopt.getOptArg();
					break;
					}
				case 'm':manifestFile=new File(getopt.getOptArg());break;
				default: 
					{
					switch(handleOtherOptions(c, getopt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default: break;
						}
					}
				}
			}
		if(record_per_file<0L)
			{
			System.err.println("Undefined number of record per file");
			return -1;
			}
		
		SAMFileWriter sfw=null;
		SAMRecordIterator iter=null;
		SAMFileReader samFileReader=null;
		PrintWriter manifest=new PrintWriter(new NullOuputStream());
		try
			{
			if(getopt.getOptInd()==args.length)
				{
				info("reading from stdin.");
				samFileReader=new SAMFileReader(System.in);
				}
			else if(getopt.getOptInd()+1==args.length)
				{
				File file=new File(args[getopt.getOptInd()]);
				info("REading "+file);
				samFileReader=new SAMFileReader(file);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			samFileReader.setValidationStringency(ValidationStringency.SILENT);
			SAMFileHeader header=samFileReader.getFileHeader();
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setCommandLine(this.getProgramCommandLine());
			prg.setProgramName(this.getProgramName());
			prg.setProgramVersion(this.getVersion());
			
			int split_file_number=0;
			long nReads=0L;
			iter=samFileReader.iterator();
			
			if(manifestFile!=null)
				{
				manifest=new PrintWriter(manifestFile);
				}
			
			while(iter.hasNext())
				{
				SAMRecord rec=iter.next();
				++nReads;
				if(sfw==null)
					{
					split_file_number++;
					String pathname=this.prefix+"."+String.format("%0"+suffix_length+"d", split_file_number)+".bam";
					File out=new File(pathname);
					info("Opening "+out);
					manifest.write(pathname);
					manifest.write("\t"+(nReads)+"\t");
					
					SAMFileHeader header2=header.clone();
					header2.addComment("SPLIT:"+split_file_number);
					header2.addComment("SPLIT:Starting from Read"+nReads);
					
					sfw=swfactory.make(header2, out);
					}
				sfw.addAlignment(rec);
				
				if(nReads%record_per_file==0)
					{
					info("Closing "+sfw);
					sfw.close();
					manifest.write((nReads)+"\n");
					sfw=null;
					}
				
				}
			if(sfw!=null)
				{
				sfw.close();
				manifest.write((nReads)+"\n");
				}
			manifest.flush();
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(sfw);
			CloserUtil.close(iter);
			CloserUtil.close(samFileReader);
			}
		return 0;
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Biostar90204().instanceMainWithExit(args);

	}

}
