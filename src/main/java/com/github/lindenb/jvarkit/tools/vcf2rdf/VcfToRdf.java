package com.github.lindenb.jvarkit.tools.vcf2rdf;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;



import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.rdf.RDFVcfWriter;

public class VcfToRdf extends AbstractVCFFilter2
	{
	private VcfToRdf()
		{
		}
	
	@Override
	public String getProgramDescription() {
		return "convert VCF to RDF";
		}
	
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfToRdf";
		}
	
	@Override
	public void printOptions(java.io.PrintStream out)
		{
		super.printOptions(out);
		}
	

	@Override
	protected VariantContextWriter createVariantContextWriter(File OUT) throws IOException
		{
		if(OUT==null)
			{
			info("writing to stdout");
			return new RDFVcfWriter(System.out);
			}
		else if(OUT.getName().endsWith(".gz"))
			{
			info("writing to "+OUT+" + gzip");
			GZIPOutputStream zout=new GZIPOutputStream(new FileOutputStream(OUT));
			return new RDFVcfWriter(zout);
			}
		else
			{
			info("writing to "+OUT);
			return new RDFVcfWriter(new FileOutputStream(OUT));
			}
		}

	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
		throws IOException
		{
		out.writeHeader(in.getHeader());
		while(in.hasNext())
			{
			out.add(in.next());
			}
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
					switch(handleOtherOptions(c, opt, null))
						{
						case EXIT_FAILURE: return -1;
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
		new VcfToRdf().instanceMainWithExit(args);
		}
	}
