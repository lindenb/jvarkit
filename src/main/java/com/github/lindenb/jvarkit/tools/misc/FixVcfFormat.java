package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.regex.Pattern;

import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class FixVcfFormat extends AbstractCommandLineProgram
	{
	private FixVcfFormat()
		{
		
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/FixVcfFormat";
		}
	
	@Override
	public String getProgramDescription() {
		return "Fix PL format in VCF. Problem is described in http://gatkforums.broadinstitute.org/discussion/3453";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+""))!=-1)
			{
			switch(c)
				{
				default:
					{
					switch(handleOtherOptions(c, opt, args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		long n_var=0L;
		BufferedReader in=null;
		long n_fix=0L;
		int report_mismatch_sample_call=0;
		try
			{
			Pattern tab=Pattern.compile("[\t]");
			Pattern colon=Pattern.compile("[\\:]");
			Pattern comma=Pattern.compile("[,]");
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				in=new BufferedReader(new InputStreamReader(System.in));
				}
			else if(opt.getOptInd()+1==args.length)
				{
				
				String filename=args[opt.getOptInd()];
				info("Reading from "+filename);
				in=IOUtils.openURIForBufferedReading(filename);
				}
			else
				{
				error("Illegal number of arguments.");
				return -1;
				}
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.startsWith("#"))
					{
					if(line.startsWith("#CHROM"))
						{
						System.out.println("##FixVcfFormatCmd="+getProgramCommandLine().replace('\n', ' '));
						System.out.println("##FixVcfFormatVersion="+getVersion());
						}
					System.out.println(line);
					continue;
					}
				if(++n_var%10000==0)
					{
					info("Variant:"+n_var+" fix:"+n_fix);
					}
				String tokens[]=tab.split(line);
				if(tokens.length<9)
					{
					warning("not enought column in "+line);
					System.out.println(line);
					continue;
					}
				String formats[]=colon.split(tokens[8]);
				int PL_index=-1;
				for(int i=0;i< formats.length;++i)
					{
					if(formats[i].equals("PL"))
						{
						PL_index=i;
						break;
						}
					}
				for(int i=0;i< 9;++i)
					{
					if(i>0) System.out.print("\t");
					System.out.print(tokens[i]);
					}
				for(int sample=9;sample< tokens.length;++sample)
					{
					System.out.print("\t");
					if(tokens[sample].equals("."))
						{
						System.out.print(tokens[sample]);
						continue;
						}
					
					String calls[]=colon.split(tokens[sample]);
					if(calls.length>formats.length)
						{
						error("not same number of columns between FORMAT and call:"+tokens[8]+" vs "+tokens[sample]);
						return -1;
						}
					else if(calls.length<formats.length)
						{
						if(report_mismatch_sample_call<10)
							{
							warning("not same number of columns between FORMAT and call:"+tokens[8]+" vs "+tokens[sample]);
							}
						report_mismatch_sample_call++;
						}
					
					for(int i=0;i< calls.length;++i)
						{
						if(i>0) System.out.print(':');
						if(i==PL_index && !calls[i].equals("."))
							{
							String pl_values[]=comma.split(calls[i]);
							for(int j=0;j< pl_values.length;++j)
								{
								if(j>0) System.out.print(",");
								if(pl_values[j].equals("."))
									{
									System.out.print("0");
									++n_fix;
									}
								else
									{
									System.out.print(pl_values[j]);
									}
								}
							}
						else
							{
							System.out.print(calls[i]);
							}
						}
					//FORMAT missing
					for(int i=calls.length;i< formats.length ; ++i)
						{
						if(i>0) System.out.print(':');
						System.out.print(".");
						}
					
					}
				System.out.println();
				if(System.out.checkError()) break;
				}
			info("Number of FIX: "+n_fix);
			return 0;
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
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FixVcfFormat().instanceMainWithExit(args);

	}

}
