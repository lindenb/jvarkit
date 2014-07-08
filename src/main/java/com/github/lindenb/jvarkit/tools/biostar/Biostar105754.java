package com.github.lindenb.jvarkit.tools.biostar;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.regex.Pattern;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

public class Biostar105754 extends AbstractCommandLineProgram
	{
	private PrintWriter out;
	private org.broad.igv.bbfile.BBFileReader bbFileReader;
	private final long EXTEND_SHIFT=1000000;//
	private final long MAX_CHROM_END=Integer.MAX_VALUE-EXTEND_SHIFT;
	private Biostar105754()
		{
		}
	@Override
	public String getProgramDescription() {
		return "bigwig : peak distance from specific genomic region";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/Biostar105754";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -B (bigwig file) Required.");
		super.printOptions(out);
		}
	
	private static int distance(int start1,int end1, int start2,int end2)
		{
		if(end1< start2)
			{
			return start2-end1;
			}
		else if(end2< start1)
			{
			return start1-end2;
			}
		int d=Math.abs(start1-start2);
		d=Math.min(d,Math.abs(start1-end2));
		d=Math.min(d,Math.abs(end1-start2));
		d=Math.min(d,Math.abs(end1-end2));
		return d;
		}
	
	private void run(LineIterator r)
		throws IOException
		{
		Pattern tab=Pattern.compile("[\t]");
		while(r.hasNext())
			{
			String line=r.next();
			if(line.startsWith("#"))
				{
				continue;
				}
			String tokens[]=tab.split(line,4);
			if(tokens.length<3)
				{
				System.err.println("Bad BED line: "+line);
				continue;
				}
			String chrom=tokens[0];
			int chromStart0=Integer.parseInt(tokens[1]);
			int chromEnd0=Integer.parseInt(tokens[2]);
			if(chrom.isEmpty() || chromStart0<0L || chromEnd0<chromStart0)
				{
				System.err.println("Bad BED line: "+line);
				continue;
				}
			
			//extends bed area until something was found
			int chromStart=chromStart0;
			int chromEnd=chromEnd0;
			
			for(;;)
				{
				BigWigIterator iter=this.bbFileReader.getBigWigIterator(
						chrom,
						chromStart,
						chrom,
						chromEnd,
						false);
				if(iter!=null)
					{
					WigItem best=null;
					while(iter.hasNext())
						{
						WigItem wigItem=iter.next();
						if(best==null || 
								distance(chromStart,chromEnd,best.getStartBase(),best.getEndBase()) >
								distance(chromStart,chromEnd,wigItem.getStartBase(),wigItem.getEndBase())
								)
							{
							best=wigItem;
							}
						}
					if(best!=null)
						{
						this.out.print(best.getChromosome());
						this.out.print("\t");
						this.out.print(best.getStartBase());
						this.out.print("\t");
						this.out.print(best.getEndBase());
						this.out.print("\t");
						this.out.print(best.getWigValue());
						this.out.print("\t");
						this.out.print(line);
						this.out.println();
						break;
						}
					}
				//extend bed area
				long start2=chromStart-EXTEND_SHIFT;
				long end2=chromEnd+EXTEND_SHIFT;
				if(start2<0) start2=0;
				if(end2>MAX_CHROM_END) end2=MAX_CHROM_END;
				//too wide, break loop
				if(start2==0 && end2==MAX_CHROM_END)
					{
					System.err.println("no data found for\t"+line);
					break;
					}
				chromStart=(int)start2;
				chromEnd=(int)end2;
				}
			}
		}
	
	@Override
	public int doWork(String[] args)
		{
		String bigWigFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"B:"))!=-1)
			{
			switch(c)
				{
				case 'B': bigWigFile= opt.getOptArg();break;
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
		if(bigWigFile==null)
			{
			error("Big wig file undefined");
			return -1;
			}
		
		try
			{
			info("Opening "+bbFileReader);
			this.bbFileReader=new BBFileReader(bigWigFile);
			if(!this.bbFileReader.isBigWigFile())
				{
				error("File "+bigWigFile+" is not a bigwig file");
				return -1;
				}
			this.out=new PrintWriter(System.out);
			
			if(opt.getOptInd()==args.length)
				{
				info("Reading BED from stdin");
				LineIterator r=IOUtils.openStdinForLineIterator();
				run(r);
				CloserUtil.close(r);
				}
			else
				{
				for(int i=opt.getOptInd();
						i< args.length && !out.checkError()
						;++i)
					{
					String filename=args[i];
					info("Reading BED from "+filename);
					LineIterator r=IOUtils.openURIForLineIterator(filename);
					run(r);
					CloserUtil.close(r);

					}
				}
			this.out.flush();
			this.out.close();
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(bbFileReader);
			CloserUtil.close(this.out);
			}
		}
	public static void main(String[] args) {
		new Biostar105754().instanceMainWithExit(args);
	}
	}
