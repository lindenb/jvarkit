/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.io.FileOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.cmdline.Option;
import htsjdk.samtools.cmdline.Usage;
import htsjdk.samtools.util.Log;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.tools.vcf2sql.VcfToSql;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
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
			
			for(SnpEffPredictionParser.SnpEffPrediction pred:parser.getPredictions(ctx))
				{
				w.writeStartElement(getKey());
				
				simpleTag2(w,"altAA",pred.getAltAminoAcid());
				simpleTag2(w,"posAA",pred.getAminoAcidPosition());
				simpleTag2(w,"ensGene",pred.getEnsemblGene());
				simpleTag2(w,"ensProtein",pred.getEnsemblProtein());
				simpleTag2(w,"ensTranscript",pred.getEnsemblTranscript());
				simpleTag2(w,"gene",pred.getGeneName());
				simpleTag2(w,"refAA",pred.getReferenceAminoAcid());
				for(SequenceOntologyTree.Term t:pred.getSOTerms())
					{
					simpleTag2(w,"so",t.getAcn());
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
		private VepPredictionParser parser=null;
		@Override
		public void handle(VCFHeader header, XMLStreamWriter w,
				VariantContext ctx) throws XMLStreamException
			
			{
			if(parser==null) parser=new VepPredictionParser(header);
			
			for(VepPredictionParser.VepPrediction pred:parser.getPredictions(ctx))
				{
				w.writeStartElement(getKey());
				
				simpleTag2(w,"altAA",pred.getAltAminoAcid());
				simpleTag2(w,"posAA",pred.getAminoAcidPosition());
				simpleTag2(w,"ensGene",pred.getEnsemblGene());
				simpleTag2(w,"ensProtein",pred.getEnsemblProtein());
				simpleTag2(w,"ensTranscript",pred.getEnsemblTranscript());
				simpleTag2(w,"gene",pred.getGeneName());
				simpleTag2(w,"exon",pred.getExon());
				simpleTag2(w,"hgnc",pred.getHGNC());
				simpleTag2(w,"refAA",pred.getReferenceAminoAcid());
				for(SequenceOntologyTree.Term t:pred.getSOTerms())
					{
					simpleTag2(w,"so",t.getAcn());
					}
				
				w.writeEndElement();
				}
			}
		
		@Override
		public String getKey() {
			return VepPredictionParser.getDefaultTag();
			}
		}
	
	private static void simpleTag2(XMLStreamWriter w,String tag,Object content) throws XMLStreamException
		{
		if(content==null || content.toString().isEmpty()) return;
		simpleTag1(w,tag,content);
		}
	private static void simpleTag1(XMLStreamWriter w,String tag,Object content) throws XMLStreamException
		{
		
		w.writeStartElement(tag);
		w.writeCharacters(String.valueOf(content));
		w.writeEndElement();
		}
	
	
}
