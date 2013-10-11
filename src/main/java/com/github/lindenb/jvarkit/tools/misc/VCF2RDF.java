package com.github.lindenb.jvarkit.tools.misc;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;

import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.rdf.RDFVcfWriter;

public class VCF2RDF extends AbstractVCFFilter
	{
	private static final Log LOG=Log.getInstance(VCF2RDF.class);
	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Convert vcf to RDF. ";
	
	public VCF2RDF()
		{
		}
	@Override
	protected VariantContextWriter createVariantContextWriter() throws IOException
		{
		if(OUT==null)
			{
			LOG.info("writing to stdout");
			return new RDFVcfWriter(System.out);
			}
		else if(OUT.getName().endsWith(".gz"))
			{
			LOG.info("writing to "+OUT+" + gzip");
			GZIPOutputStream zout=new GZIPOutputStream(new FileOutputStream(OUT));
			return new RDFVcfWriter(zout);
			}
		else
			{
			LOG.info("writing to "+OUT);
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
	
	public static void main(String[] args)
		{
		new VCF2RDF().instanceMainWithExit(args);
		}

}
