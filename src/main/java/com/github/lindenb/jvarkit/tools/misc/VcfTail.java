package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedList;


import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfTail extends AbstractVCFFilter2
	{
	private long count=10L;
	private VcfTail()
		{
		}
	

	@Override
	public String getProgramDescription() {
		return "Print last lines of a VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfTail";
		}
	
	

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		VCFHeader header=in.getHeader();
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
		out.writeHeader(header);
		LinkedList<VariantContext> L=new LinkedList<VariantContext>();
		while(in.hasNext() && L.size()< this.count)
			{	
			L.add(in.next());
			}
		while(in.hasNext())
			{
			L.add(in.next());
			L.removeFirst();
			}
		for(VariantContext ctx:L)
			{
			out.add(ctx);
			}
		}
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -n (int) output size. Optional. Default:"+this.count);
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "n:"))!=-1)
			{
			switch(c)
				{
				case 'n': this.count=Math.max(0, Long.parseLong(opt.getOptArg())); break;
				default: 
					{
					switch(handleOtherOptions(c, opt))
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
		new VcfTail().instanceMainWithExit(args);
		}
	}
