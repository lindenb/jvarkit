package com.github.lindenb.jvarkit.tools.vcf2rdf;

import java.io.File;
import java.io.IOException;
import java.net.URI;



import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.rdf.RDFVcfWriter;

public class VcfToRdf extends AbstractCommandLineProgram
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
	


	protected void doWork(VcfIterator in, RDFVcfWriter out,URI src)
		throws IOException
		{
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader().getSequenceDictionary());
		out.writeHeader(in.getHeader(),src);
		while(in.hasNext())
			{
			VariantContext ctx=in.next();
			progress.watch(ctx.getChr(), ctx.getStart());
			out.add(ctx);
			}
		progress.finish();
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
		VcfIterator in=null;
		RDFVcfWriter w=null;
		try
			{
		
			if(opt.getOptInd()==args.length)
				{
				w=new RDFVcfWriter(System.out);
				in= VCFUtils.createVcfIteratorStdin();
				doWork(in,w,null);
				in.close();
				w.close();
				}
			else
				{
				w=new RDFVcfWriter(System.out);
				w.writeStartDocument();
				for(int i=opt.getOptInd();i< args.length;++i )
					{
					URI source=null;
					in=VCFUtils.createVcfIterator(args[i]);
					if(IOUtils.isRemoteURI(args[i]))
						{
						source=URI.create(args[i]);
						}
					else
						{
						source=new File(args[i]).toURI();
						}
					doWork(in,w,source);
					in.close();
					}
				w.close();
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
			CloserUtil.close(in);
			CloserUtil.close(w);
			}
		
		}
	public static void main(String[] args)
		{
		new VcfToRdf().instanceMainWithExit(args);
		}
	}
