package com.github.lindenb.jvarkit.tools.misc;
/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.




*/
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Collection;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

public class FixVcfFormat extends AbstractFixVcfFormat
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FixVcfFormat.class);

	
	
	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFixVcfFormat.AbstractFixVcfFormatCommand
	 	{		
		@SuppressWarnings("resource")
		@Override
		protected Collection<Throwable> call(String inputName) throws Exception
				{
				long n_var=0L;
				BufferedReader in=null;
				long n_fix=0L;
				PrintStream out = null;
				int report_mismatch_sample_call=0;
				try
					{
					Pattern tab=Pattern.compile("[\t]");
					Pattern colon=Pattern.compile("[\\:]");
					Pattern comma=Pattern.compile("[,]");
					if(inputName==null)
						{
						LOG.info("Reading from stdin");
						in=new BufferedReader(new InputStreamReader(System.in));
						}
					else
						{
						LOG.info("Reading from "+inputName);
						in=IOUtils.openURIForBufferedReading(inputName);
						}
					out  = openFileOrStdoutAsPrintStream();
					String line;
					while((line=in.readLine())!=null)
						{
						if(line.startsWith("#"))
							{
							if(line.startsWith("#CHROM"))
								{
								out.println("##FixVcfFormatCmd="+getProgramCommandLine().replace('\n', ' '));
								out.println("##FixVcfFormatVersion="+getVersion());
								}
							out.println(line);
							continue;
							}
						if(++n_var%10000==0)
							{
							LOG.info("Variant:"+n_var+" fix:"+n_fix);
							}
						String tokens[]=tab.split(line);
						if(tokens.length<9)
							{
							LOG.warn("not enough column in "+line);
							out.println(line);
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
							if(i>0) out.print("\t");
							out.print(tokens[i]);
							}
						for(int sample=9;sample< tokens.length;++sample)
							{
							out.print("\t");
							if(tokens[sample].equals("."))
								{
								out.print(tokens[sample]);
								continue;
								}
							
							String calls[]=colon.split(tokens[sample]);
							if(calls.length>formats.length)
								{
								return wrapException("not same number of columns between FORMAT and call:"+tokens[8]+" vs "+tokens[sample]);
								}
							else if(calls.length<formats.length)
								{
								if(report_mismatch_sample_call<10)
									{
									LOG.warn("not same number of columns between FORMAT and call:"+tokens[8]+" vs "+tokens[sample]);
									}
								report_mismatch_sample_call++;
								}
							
							for(int i=0;i< calls.length;++i)
								{
								if(i>0) out.print(':');
								if(i==PL_index && !calls[i].equals("."))
									{
									String pl_values[]=comma.split(calls[i]);
									for(int j=0;j< pl_values.length;++j)
										{
										if(j>0) out.print(",");
										if(pl_values[j].equals("."))
											{
											out.print("0");
											++n_fix;
											}
										else
											{
											out.print(pl_values[j]);
											}
										}
									}
								else
									{
									out.print(calls[i]);
									}
								}
							//FORMAT missing
							for(int i=calls.length;i< formats.length ; ++i)
								{
								if(i>0) out.print(':');
								out.print(".");
								}
							
							}
						out.println();
						if(out.checkError()) break;
						}
					out.flush();
					LOG.info("Number of FIX: "+n_fix);
					return RETURN_OK;
					}
				catch(Exception err)
					{
					LOG.error(err);
					return wrapException(err);
					}
				finally
					{
					CloserUtil.close(out);
					CloserUtil.close(in);
					}
				}
	 	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FixVcfFormat().instanceMainWithExit(args);
	}

}
