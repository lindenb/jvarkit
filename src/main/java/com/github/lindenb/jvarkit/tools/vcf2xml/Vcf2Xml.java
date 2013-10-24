/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.GZIPOutputStream;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.util.Log;

import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.tools.vcf2sql.VcfToSql;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.XMLVcfWriter;
import com.github.lindenb.jvarkit.util.vcf.predictions.SnpEffPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;

/**
 * @author lindenb
 *
 */
public class Vcf2Xml extends AbstractVCFFilter
	{
    private static Log LOG=Log.getInstance(VcfToSql.class); 

	@Usage(programVersion="1.0")
	public String USAGE=getStandardUsagePreamble()+" Convert VCF to XML";

    @Option(shortName="VEP",doc="Use  and explode VEP predictions",optional=true)
	public boolean USE_VEP=true;
    @Option(shortName="SNPEFF",doc="Use and explode SNPEFF predictions",optional=true)
	public boolean USE_SNPEFF=true;
	
	/* (non-Javadoc)
	 * @see com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter#doWork(com.github.lindenb.jvarkit.util.vcf.VcfIterator, org.broadinstitute.variant.variantcontext.writer.VariantContextWriter)
	 */
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException {
			out.writeHeader(in.getHeader());
			while(in.hasNext())
				{
				out.add(in.next());
				}
			
			}

	protected VariantContextWriter createVariantContextWriter() throws IOException
		{
		XMLVcfWriter vc=null;
		if(OUT==null)
			{
			LOG.info("writing to stdout");
			vc= new XMLVcfWriter(System.out);
			}
		else if(OUT.getName().endsWith(".gz"))
			{
			LOG.info("writing to "+OUT+" as gz file.");
			GZIPOutputStream bcos=new GZIPOutputStream(new FileOutputStream(OUT));
			vc=  new XMLVcfWriter(bcos);
			}
		else
			{
			LOG.info("writing to "+OUT);
			vc= new XMLVcfWriter(new FileOutputStream(OUT));
			}
		if(USE_SNPEFF)
			{
			vc.putInfoHandler(new SnpEffInfoHandler());
			}
		if(USE_VEP)
			{
			vc.putInfoHandler(new VepInfoHandler());
			}
		return vc;
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Vcf2Xml().instanceMainWithExit(args);
		}

	private static class SnpEffInfoHandler implements XMLVcfWriter.XMLInfoHandler
		{	
		private SnpEffPredictionParser parser=null;
		@Override
		public void handle(VCFHeader header, XMLStreamWriter w,
				VariantContext ctx) throws XMLStreamException {
			if(parser==null) parser=new SnpEffPredictionParser(header);
			
			for(Map<String,String> pred:parser.split(ctx))
				{
				w.writeStartElement(getKey());
				for(String k:pred.keySet())
					{
					w.writeStartElement(k);
					w.writeCharacters(pred.get(k));
					w.writeEndElement();
					}
				
				w.writeEndElement();
				}
			}
		@Override
		public String getKey() {
			return SnpEffPredictionParser.getDefaultTag();
			}
		}
	private static class VepInfoHandler  implements XMLVcfWriter.XMLInfoHandler
		{	
		private Pattern amp=Pattern.compile("[&]");
		private VepPredictionParser parser=null;
		@Override
		public void handle(VCFHeader header, XMLStreamWriter w,
				VariantContext ctx) throws XMLStreamException
			
			{
			if(parser==null) parser=new VepPredictionParser(header);
			
			for(Map<String,String> pred:parser.split(ctx))
				{
				w.writeStartElement(getKey());
				for(String k:pred.keySet())
					{
					String array[]=new String[]{pred.get(k)};
					if(k.equals("Consequence"))
						{
						array=amp.split(array[0]);
						}
					
					for(String item:array)
						{
						if(item.isEmpty()) continue;
						w.writeStartElement(k);
						w.writeCharacters(item);
						w.writeEndElement();
						}
					}
				
				w.writeEndElement();
				}
			}
		
		@Override
		public String getKey() {
			return VepPredictionParser.getDefaultTag();
			}
		}

}
