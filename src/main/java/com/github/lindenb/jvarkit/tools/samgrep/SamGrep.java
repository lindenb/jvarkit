package com.github.lindenb.jvarkit.tools.samgrep;


import java.io.BufferedReader;
import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SamWriterFactory;


import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloserUtil;

public class SamGrep extends AbstractCommandLineProgram
	{
    private SamGrep()
    	{
    	}
    
	@Override
	public String getProgramDescription()
		{
		return "grep read-names in a bam file";
		}
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/SamGrep";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		
		out.println(" -f (file) file containing a list of read names..");
		out.println(" -R (name) add the read.");
		out.println(" -n (int) when found, remove the read from the list of names when found more that 'n' time (increase speed)");
		out.println(" -V  invert");

		out.println(" -o (filename) output file. default: stdout.");
		out.println(" -X if -o used, continue to output original input to stdout.");
		out.println(" -c (int) compression level");
		out.println(" -b force binary");
		super.printOptions(out);
		}

    
	
	@Override
	public int doWork(String[] args)
		{
		boolean divertToStdout=false;
		int n_before_remove=-1;
		Map<String,Integer> readNames=new HashMap<String,Integer>(); 
		SamWriterFactory swf=SamWriterFactory.newInstance();
		File fileout=null;
		boolean inverse=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:bc:f:R:n:X"))!=-1)
			{
			switch(c)
				{
				case 'X': divertToStdout=true;break;
				case 'n': n_before_remove=Integer.parseInt(opt.getOptArg()); break;
				case 'V': inverse=true;break;
				case 'R': readNames.put(opt.getOptArg(),0);break;
				case 'f':
					{
					BufferedReader in=null;
					try
						{
						in=IOUtils.openURIForBufferedReading(opt.getOptArg());
				    	String line;
				    	while((line=in.readLine())!=null)
				    		{
				    		line=line.trim();
				    		if(line.isEmpty()) continue;
				    		readNames.put(line,0);
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
				case 'o': fileout=new File(opt.getOptArg());break;
				case 'b': swf.setBinary(true);break;
				case 'c': swf.setCompressionLevel(Integer.parseInt(opt.getOptArg()));break;
				default:
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(readNames.isEmpty())
    		{
    		warning("no read found.");
    		}
		
		SAMFileWriter sfw=null;
		SAMFileWriter samStdout=null;
		SAMFileReader sfr=null;
		try
			{
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading sfomr stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				info("Reading from "+filename);
				sfr=new SAMFileReader(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.LENIENT);
			SAMFileHeader header=sfr.getFileHeader().clone();
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setProgramName(getProgramName());
			prg.setProgramVersion(getVersion());
			prg.setCommandLine(getProgramCommandLine());
			
			
			if(fileout==null)
				{
				sfw=swf.make(header);
				}
			else
				{
				if(divertToStdout) samStdout=swf.make(sfr.getFileHeader());
				sfw=swf.make(header,fileout);
				}
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				boolean keep=false;
				SAMRecord rec=iter.next();
				if(samStdout!=null) samStdout.addAlignment(rec);
				Integer count=readNames.get(rec.getReadName());
				if(count!=null)
					{
					keep=true;
					}
				if(inverse) keep=!keep;
				if(keep)
					{
					sfw.addAlignment(rec);
					}
				
				if(n_before_remove!=-1 && !inverse && keep)
					{
					count++;
					if(count>=n_before_remove)
						{
						readNames.remove(rec.getReadName());
						if(samStdout==null && readNames.isEmpty()) break;
						}
					else
						{
						readNames.put(rec.getReadName(),count);
						}
					}
				
				}
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(samStdout);
			CloserUtil.close(sfw);
			CloserUtil.close(sfr);
			}
		return 0;
		}

    public static void main(final String[] argv)
		{
	    new SamGrep().instanceMainWithExit(argv);
		}	

	}
