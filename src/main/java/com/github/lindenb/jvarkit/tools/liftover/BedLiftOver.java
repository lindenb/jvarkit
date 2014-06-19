package com.github.lindenb.jvarkit.tools.liftover;

import java.io.File;
import java.io.PrintStream;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.CloserUtil;



public class BedLiftOver extends AbstractCommandLineProgram
	{
	private LiftOver liftOver=null;
	@Override
	public String getProgramDescription() {
		return "Lift-over a BED file.";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/BedLiftOver";
		}

	
	private void scan(LineIterator r)
		{
		Pattern tab=Pattern.compile("\t");
		while(r.hasNext())
			{
			String line=r.next();
			if(line.startsWith("#") || line.trim().isEmpty()) continue;
			String tokens[]=tab.split(line, 4);
			if(tokens.length<3) continue;
			String chrom=tokens[0];
			int start0=Integer.parseInt(tokens[1]);
			int end0=Integer.parseInt(tokens[2]);
			Interval srcInterval=new Interval(chrom, start0+1,end0);
			Interval dest=this.liftOver.liftOver(srcInterval);
			if(dest!=null)
				{
				System.out.print(dest.getSequence());
				System.out.print('\t');
				System.out.print(dest.getStart()-1);
				System.out.print('\t');
				System.out.print(dest.getEnd());
				System.out.print('\t');
				System.out.print(line);
				System.out.println();
				}
			else
				{
				List<LiftOver.PartialLiftover> L=this.liftOver.diagnosticLiftover(srcInterval);
				if(L!=null && !L.isEmpty())
					{	
					for(LiftOver.PartialLiftover plo:L)
						{	
						System.out.print("#");
						System.out.print(tokens[0]);
						System.out.print('\t');
						System.out.print(tokens[1]);
						System.out.print('\t');
						System.out.print(tokens[2]);
						System.out.print('\t');
						System.out.println(plo);
						}
					}
				else
					{
					System.out.print("#");
					System.out.print(line);
					}
				}			
			}
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -f (chain-file) LiftOver file. Required.");
		out.println(" -m (double) lift over min-match. default:"+LiftOver.DEFAULT_LIFTOVER_MINMATCH);
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		double minMatch=LiftOver.DEFAULT_LIFTOVER_MINMATCH;
		File liftOverFile=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "f:m:"))!=-1)
			{
			switch(c)
				{
				case 'f': liftOverFile=new File(opt.getOptArg()); break;
				case 'm': minMatch=Double.parseDouble(opt.getOptArg()); break;
				default: 
					{
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE:return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(liftOverFile==null)
			{
			error("LiftOver file is undefined.");
			return -1;
			}
		this.liftOver=new LiftOver(liftOverFile);
		this.liftOver.setLiftOverMinMatch(minMatch);
		try
			{
			if(opt.getOptInd()==args.length)
				{
				LineIterator r=IOUtils.openStdinForLineIterator();
				scan(r);
				CloserUtil.close(r);
				}
			else
				{
				for(int optind=opt.getOptInd();optind< args.length;++optind)
					{
					String filename=args[optind];
					LineIterator r=IOUtils.openURIForLineIterator(filename);
					scan(r);
					CloserUtil.close(r);
					}
				}
			return 0;
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		
		}

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new BedLiftOver().instanceMainWithExit(args);
		}

	}
