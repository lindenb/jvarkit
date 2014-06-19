package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class DownSampleVcf extends AbstractVCFFilter2
	{
	private int reservoir_size=10;
	private long seed=System.currentTimeMillis();
	private DownSampleVcf()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "DownSample a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/DownSampleVcf";
		}
	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		Random rand=new Random(this.seed);
		List<VariantContext>  buffer=new ArrayList<VariantContext>(this.reservoir_size);
		VCFHeader header=in.getHeader();
		header.addMetaDataLine(new VCFHeaderLine("DownSampleVcf_CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine("DownSampleVcf_Version",String.valueOf(getVersion())));
		out.writeHeader(header);
		if(this.reservoir_size!=0)
			{
			while(in.hasNext())
				{	
				if(buffer.size() < this.reservoir_size)
					{
					buffer.add(in.next());
					}
				else
					{
					buffer.set(rand.nextInt(buffer.size()), in.next());
					}
				}
			}
		for(VariantContext ctx:buffer)
			{
			out.add(ctx);
			}
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -N (long) random seed. Optional.");
		out.println(" -n (int) output size. Optional");
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "N:n:"))!=-1)
			{
			switch(c)
				{
				case 'N': seed=Long.parseLong(opt.getOptArg()); break;
				case 'n': reservoir_size=Math.max(0, Integer.parseInt(opt.getOptArg())); break;
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
		
		return doWork(opt.getOptInd(), args);
		}

	public static void main(String[] args)
		{
		new DownSampleVcf().instanceMainWithExit(args);
		}
	}
