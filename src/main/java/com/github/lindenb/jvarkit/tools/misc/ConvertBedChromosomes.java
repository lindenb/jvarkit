package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReaderUtil;

import net.sf.samtools.util.CloserUtil;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class ConvertBedChromosomes
	extends AbstractCommandLineProgram
	{
	private Map<String,String> customMapping=new HashMap<String,String>();
	private Set<String> unmappedChromosomes=new HashSet<String>();
	private int chromColumn0=0;
	private ConvertBedChromosomes()
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
			return null;
			}
		return newname;
		}
	
	@SuppressWarnings("resource")
	protected int doWork(InputStream in,PrintStream out)
			throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		LineIterator lr=new LineIteratorImpl(LineReaderUtil.fromBufferedStream(in));
		while(lr.hasNext())
			{	
			String line=lr.next();
			String tokens[]=tab.split(line, (chromColumn0+2));
			if(chromColumn0 >=tokens.length) throw new IOException("Bad BED line : "+line+" extected at least "+(chromColumn0+2)+" columns");
			String chrom=convertName(tokens[chromColumn0]);
			if(chrom==null) continue;
			for(int i=0;i< tokens.length;++i)
				{
				if(i>0) out.print("\t");
				out.print(i==chromColumn0?chrom:tokens[i]);
				}
			out.println();
			}
		
		return 0;
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BedRenameChromosomes";
		}
	
	@Override
	public String getProgramDescription() {
		return "Convert the names of the chromosomes in a Bed file.";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -f (file) load a custom name mapping. Format (chrom-source\\tchrom-dest\\n)+");
		out.println(" -c (int) 1-based chromosome column: default: 1");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"f:c:"))!=-1)
			{
			switch(c)
				{
				case 'c':
						this.chromColumn0=Integer.parseInt(opt.getOptArg())-1;
						if(this.chromColumn0<0)
							{
							error("bad chromosome index (<1): "+opt.getOptArg());
							return -1;
							}
						break;
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
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		
		if(customMapping.isEmpty())
			{
			error("No custom mapping defined");
			}
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("reading stdin");
				doWork(System.in, System.out);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
				
					{
					info("opening "+args[i]);
					InputStream in=IOUtils.openURIForReading(args[i]);
					doWork(in, System.out);
					CloserUtil.close(in);
					}
				}
			if(!unmappedChromosomes.isEmpty())
				{
				warning("Unmapped chromosomes:"+unmappedChromosomes);
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		
		}
	

	public static void main(String[] args)
		{
		new ConvertBedChromosomes().instanceMainWithExit(args);
		}
	}
