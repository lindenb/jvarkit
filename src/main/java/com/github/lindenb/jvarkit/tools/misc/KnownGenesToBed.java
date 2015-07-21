package com.github.lindenb.jvarkit.tools.misc;

import java.io.PrintStream;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

public class KnownGenesToBed extends AbstractCommandLineProgram
	{
	private PrintStream out=System.out;
	private boolean print_introns=true;
	private boolean print_utrs=true;
	private boolean print_cds=true;
	private boolean print_exons=true;
	private boolean print_transcripts=true;
	
	@Override
	public String getProgramDescription() {
		return "KnownGene to BED.";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/KnownGenesToBed";
		}
	
	private void print(KnownGene kg,int start,int end,String type,String name)
		{
		if(start>=end) return;
		out.print(kg.getChromosome());
		out.print('\t');
		out.print(start);
		out.print('\t');
		out.print(end);
		out.print('\t');
		out.print(kg.isPositiveStrand()?'+':'-');
		out.print('\t');
		out.print(kg.getName());
		out.print('\t');
		out.print(type);
		out.print('\t');
		out.print(name);
		out.println();
		}
	
	private void scan(LineIterator r)
		{
		Pattern tab=Pattern.compile("[\t]");
		while(r.hasNext())
			{
			if(out.checkError()) break;
			String line=r.next();
			String tokens[]=tab.split(line);
			KnownGene kg=new KnownGene(tokens);
			if(print_transcripts) print(kg,kg.getTxStart(),kg.getTxEnd(),"TRANSCRIPT",kg.getName());
			for(int i=0;i< kg.getExonCount();++i)
				{
				KnownGene.Exon exon=kg.getExon(i);
				if(print_exons) print(kg,exon.getStart(),exon.getEnd(),"EXON",exon.getName());
				
				if(print_utrs && kg.getCdsStart()>exon.getStart())
					{
					print(kg,exon.getStart(),
							Math.min(kg.getCdsStart(),exon.getEnd()),"UTR","UTR"+(kg.isPositiveStrand()?"5":"3"));
					}
				
				if(print_cds && !(kg.getCdsStart()>=exon.getEnd() || kg.getCdsEnd()<exon.getStart()))
					{
					print(kg,
							Math.max(kg.getCdsStart(),exon.getStart()),
							Math.min(kg.getCdsEnd(),exon.getEnd()),
							"CDS",exon.getName()
							);
					}
				
				KnownGene.Intron intron=exon.getNextIntron();
				if(print_introns && intron!=null)
					{
					print(kg,intron.getStart(),intron.getEnd(),"INTRON",intron.getName());
					}
				
				if(print_utrs && kg.getCdsEnd()<exon.getEnd())
					{
					print(kg,Math.max(kg.getCdsEnd(),exon.getStart()),
							exon.getEnd(),
							"UTR","UTR"+(kg.isPositiveStrand()?"3":"5"));
					}
				
				}
			}
		}
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		out.println(" -i don't print introns");
		out.println(" -u don't print utr");
		out.println(" -c don't print cds");
		out.println(" -x don't print exons");
		out.println(" -t don't print transcripts");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"iucxt"))!=-1)
			{
			switch(c)
				{
				case 'i': print_introns=false;break;
				case 'u': print_utrs=false;break;
				case 'c': print_cds=false;break;
				case 'x': print_exons=false;break;
				case 't': print_transcripts=false;break;
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
		
		LineIterator r=null;
		try
			{
			if(opt.getOptInd()==args.length)
				{
				info("Reading from stdin");
				r=IOUtils.openStdinForLineIterator();
				scan(r);
				CloserUtil.close(r);
				}
			else
				{
				for(int i=opt.getOptInd();i< args.length;++i)
					{
					String filename=args[i];
					info("Reading from "+filename);
					r=IOUtils.openURIForLineIterator(filename);
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
		finally
			{
			
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
	new KnownGenesToBed().instanceMainWithExit(args);

	}

}
