package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC
## Example
```bash
$ java -jar dist/fixvcf.jar < bad.vcf > ok.vcf
```
END_DOC

 */
@Program(
		name="fixvcfformat",
		description="Fix PL format in VCF. Problem is described in http://gatkforums.broadinstitute.org/discussion/3453")
public class FixVcfFormat extends Launcher
	{
	private static Logger LOG=Logger.build(FixVcfFormat.class).make();

	

	private FixVcfFormat()
		{
		
		}
	
	
	@Override
	public int doWork(List<String> args)
		{
		long n_var=0L;
		BufferedReader in=null;
		long n_fix=0L;
		int report_mismatch_sample_call=0;
		try
			{
			Pattern tab=Pattern.compile("[\t]");
			Pattern colon=Pattern.compile("[\\:]");
			Pattern comma=Pattern.compile("[,]");
			in  = super.openBufferedReader(oneFileOrNull(args));
			
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
					LOG.info("Variant:"+n_var+" fix:"+n_fix);
					}
				String tokens[]=tab.split(line);
				if(tokens.length<9)
					{
					LOG.warning("not enought column in "+line);
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
						LOG.error("not same number of columns between FORMAT and call:"+tokens[8]+" vs "+tokens[sample]);
						return -1;
						}
					else if(calls.length<formats.length)
						{
						if(report_mismatch_sample_call<10)
							{
							LOG.warning("not same number of columns between FORMAT and call:"+tokens[8]+" vs "+tokens[sample]);
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
			LOG.info("Number of FIX: "+n_fix);
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
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
