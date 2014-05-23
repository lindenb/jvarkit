package com.github.lindenb.jvarkit.tools.samgrep;


import java.io.BufferedReader;
import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.picard.SamFileReaderFactory;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.util.CloserUtil;

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
		out.println(" -b force binary");
		super.printOptions(out);
		}

    
	
	@Override
	public int doWork(String[] args)
		{
		boolean binary=false;
		boolean divertToStdout=false;
		int n_before_remove=-1;
		Map<String,Integer> readNames=new HashMap<String,Integer>(); 
		SAMFileWriterFactory swf=new SAMFileWriterFactory();
		File fileout=null;
		boolean inverse=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "o:bf:R:n:X"))!=-1)
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
				case 'b': binary=true;break;
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
		SamReader sfr=null;
		try
			{
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading sfomr stdin");
				sfr=SamFileReaderFactory.mewInstance().openStdin();
				}
			else if(opt.getOptInd()+1==args.length)
				{
				File filename=new File(args[opt.getOptInd()]);
				info("Reading from "+filename);
				sfr=SamFileReaderFactory.mewInstance().open(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			SAMFileHeader header=sfr.getFileHeader().clone();
			SAMProgramRecord prg=header.createProgramRecord();
			prg.setProgramName(getProgramName());
			prg.setProgramVersion(getVersion());
			prg.setCommandLine(getProgramCommandLine());
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			
			if(fileout==null)
				{
				sfw=(binary?
						swf.makeBAMWriter(header, true, System.out)
						:swf.makeSAMWriter(header, true, System.out));
				}
			else
				{
				if(divertToStdout) samStdout=(binary?
						swf.makeBAMWriter(header, true, System.out)
						:swf.makeSAMWriter(header, true, System.out));
				sfw=(binary?
						swf.makeBAMWriter(header, true,fileout)
						:swf.makeSAMWriter(header, true,fileout));
				}
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext())
				{
				boolean keep=false;
				SAMRecord rec=iter.next();
				progress.watch(rec);
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
			progress.finish();
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
