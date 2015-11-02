package com.github.lindenb.jvarkit.tools.liftover;

import java.io.PrintStream;
import java.util.Collection;
import java.util.List;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;

import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.CloserUtil;



public class BedLiftOver extends AbstractBedLiftOver
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(BedLiftOver.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBedLiftOver.AbstractBedLiftOverCommand
		{    
		private LiftOver liftOver=null;

	
	private void scan(PrintStream out,LineIterator r)
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
				out.print(dest.getContig());
				out.print('\t');
				out.print(dest.getStart()-1);
				out.print('\t');
				out.print(dest.getEnd());
				out.print('\t');
				out.print(line);
				out.println();
				}
			else
				{
				List<LiftOver.PartialLiftover> L=this.liftOver.diagnosticLiftover(srcInterval);
				if(L!=null && !L.isEmpty())
					{	
					for(LiftOver.PartialLiftover plo:L)
						{	
						out.print("#");
						out.print(tokens[0]);
						out.print('\t');
						out.print(tokens[1]);
						out.print('\t');
						out.print(tokens[2]);
						out.print('\t');
						out.println(plo);
						}
					}
				else
					{
					out.print("#");
					out.print(line);
					}
				}			
			}
		}
	
	@Override
	public Collection<Throwable> call() throws Exception {
		if(liftOverFile==null)
			{
			return wrapException("LiftOver file is undefined.");
			}
		List<String> args = getInputFiles();
		PrintStream out=null;
		try
			{
			this.liftOver=new LiftOver(liftOverFile);
			this.liftOver.setLiftOverMinMatch(minMatch);
			out=openFileOrStdoutAsPrintStream();
				
			if(args.isEmpty())
				{
				LineIterator r=IOUtils.openStreamForLineIterator(stdin());
				scan(out,r);
				CloserUtil.close(r);
				}
			else
				{
				for(String filename : args)
					{
					LOG.info("read "+filename);
					LineIterator r=IOUtils.openURIForLineIterator(filename);
					scan(out,r);
					CloserUtil.close(r);
					}
				}
			return RETURN_OK;
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			this.liftOver=null;
			CloserUtil.close(out);
			}
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
