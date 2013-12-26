package com.github.lindenb.jvarkit.tools.vcfstripannot;

import java.io.IOException;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;


import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VCFStripAnnotations extends AbstractVCFFilter2
	{
	private Set<String> KEY=new HashSet<String>();
	private boolean RESET_FILTER=false;
	
	
	@Override
	protected String getOnlineDocUrl()
		{
		return "https://github.com/lindenb/jvarkit/wiki/VCFStripAnnotations";
		}
	
	@Override
	protected String getProgramCommandLine()
		{
		return " Removes one or more field from the INFO column of a VCF.";
		}
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -k (key) remove this INFO attribute");
		out.println(" -F reset filters.");
		super.printOptions(out);
		}
	
	@Override
	protected void doWork(VcfIterator r, VariantContextWriter w)
			throws IOException
			{
			VCFHeader header=r.getHeader();
			
			
			VCFHeader h2=new VCFHeader(header.getMetaDataInInputOrder(),header.getSampleNamesInOrder());
			
			for(Iterator<VCFInfoHeaderLine> h=h2.getInfoHeaderLines().iterator();
					h.hasNext();)
				{
				VCFInfoHeaderLine vih=h.next();
				if(this.KEY.contains(vih.getID()))
					h.remove();
				}
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
			header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));
			w.writeHeader(h2);
			
			while(r.hasNext())
				{
				VariantContext ctx=r.next();
				VariantContextBuilder b=new VariantContextBuilder(ctx);
				for(String key:KEY) b.rmAttribute(key);
				if(RESET_FILTER) b.unfiltered();
				w.add(b.make());
				}		
			}
	
	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "k:F"))!=-1)
			{
			switch(c)
				{
				case 'k': this.KEY.add(opt.getOptArg()); break;
				case 'F': RESET_FILTER=true; break;
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

	
	public static void main(String[] args) throws IOException
		{
		new VCFStripAnnotations().instanceMainWithExit(args);
		}
	}
