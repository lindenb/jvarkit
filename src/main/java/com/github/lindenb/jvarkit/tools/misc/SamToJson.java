package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintWriter;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SamJsonWriter;

public class SamToJson extends AbstractCommandLineProgram
	{
	private PrintWriter out=new PrintWriter(System.out);
	@Override
	public String getProgramDescription()
		{
		return "Convert a SAM input to JSON.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println("-n add a '\\n' after each read");
		out.println("-H don't print SAM header");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		boolean print_header=true;
		boolean crlf=false;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"nH"))!=-1)
			{
			switch(c)
				{
				case 'H': print_header=false;break;
				case 'n': crlf=true; break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		SAMFileReader sfr=null;
		SamJsonWriter swf=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				sfr=new SAMFileReader(System.in);
				}
			else if(opt.getOptInd()+1==args.length)
				{	
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				sfr=new SAMFileReader(new File(filename));	
				}
			else
				{
				error(getMessageBundle("illegal.number.of.arguments"));
				return -1;
				}
			sfr.setValidationStringency(ValidationStringency.SILENT);
			swf=new SamJsonWriter(out, sfr.getFileHeader());
			swf.setAddCarriageReturn(crlf);
			swf.setPrintHeader(print_header);
			SAMRecordIterator iter=sfr.iterator();
			while(iter.hasNext() && !out.checkError())
				{
				swf.addAlignment(iter.next());
				}
			iter.close();
			swf.close();
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
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new SamToJson().instanceMainWithExit(args);

	}

}
