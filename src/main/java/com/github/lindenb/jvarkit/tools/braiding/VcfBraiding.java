package com.github.lindenb.jvarkit.tools.braiding;

import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
 * 
 * # See also:
 *    * https://visdunneright.github.io/sequence_braiding/docs/
 * 
 */
@Program(name="vcfbraiding",
	description="visualization for variants and attributes using https://visdunneright.github.io/sequence_braiding/docs/ .",
	keywords={"vcf","visualization"},
	creationDate="20201021",
	modificationDate="20201021"
	)
public class VcfBraiding extends Launcher {
	private static final Logger LOG = Logger.build(VcfBraiding.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;


	@Override
	public int doWork(final List<String> args) {
		try {
			List<String> samples = null;
			final List<VariantContext> variants = new ArrayList<VariantContext>();
			final String input = oneFileOrNull(args);

			try(VCFIterator r= super.openVCFIterator(input)) {
				final VCFHeader header = r.getHeader();
				samples = header.getGenotypeSamples();
				while(r.hasNext()) {
					final VariantContext vc = r.next();
					if(!vc.isVariant()) continue;
					variants.add(vc);
					}
				}

			final String base="https://visdunneright.github.io/sequence_braiding/docs/";
			final String svgid="vcfid";
				try(PrintStream out = super.openPathOrStdoutAsPrintStream(this.outputFile)) {
				
				XMLOutputFactory factory= XMLOutputFactory.newInstance();
				XMLStreamWriter w= factory.createXMLStreamWriter(out,"UTF-8");
				w.writeStartDocument("UTF-8","1.0");
				w.writeStartElement("html");
				w.writeAttribute("lang","en");
				w.writeCharacters("\n");
				w.writeStartElement("head");
					w.writeStartElement("title");
						w.writeCharacters(input==null?"Sequence Braiding":input);
					w.writeEndElement();//title
					
					w.writeEmptyElement("meta");
					w.writeAttribute("charset","utf-8");
	
					
					w.writeEmptyElement("meta");
					w.writeAttribute("name","description");
					w.writeAttribute("content","");
			
					w.writeEmptyElement("meta");
					w.writeAttribute("name","author");
					w.writeAttribute("content","Pierre Lindenbaum");
					
					
			
					w.writeEmptyElement("meta");
					w.writeAttribute("name","viewport");
					w.writeAttribute("content","width=device-width, initial-scale=1");
			
					w.writeStartElement("link");
					w.writeAttribute("href","//fonts.googleapis.com/css?family=Raleway:400,300,600");
					w.writeAttribute("rel","stylesheet");
					w.writeAttribute("type","text/css");
					w.writeEndElement();
			
					w.writeStartElement("link");
					w.writeAttribute("rel","stylesheet");
					w.writeAttribute("href",base + "css/normalize.css");
					w.writeEndElement();
			
					w.writeStartElement("link");
					w.writeAttribute("rel","stylesheet");
					w.writeAttribute("href",base + "css/skeleton.css");
					w.writeEndElement();
					
					w.writeStartElement("link");
					w.writeAttribute("rel","icon");
					w.writeAttribute("type","image/png");
					w.writeAttribute("href",base + "images/favicon.png");
					w.writeEndElement();
			
					w.writeStartElement("script");
					w.writeAttribute("src","https://d3js.org/d3.v5.min.js");
					w.writeEndElement();
					
				w.writeEndElement();//head
				w.writeStartElement("body");
	
				w.writeStartElement("svg");
				w.writeAttribute("id",svgid);
				w.writeAttribute("style","margin-left: 5%");
				w.writeEndElement();//svg
				w.writeEndElement();//body
				
				w.writeStartElement("script");
				w.writeAttribute("src",base+"../dist/sequence_braiding.js");
				w.writeEndElement();//script
					
				//
				w.writeStartElement("script");
				w.writeAttribute("type","text/javascript");

				w.writeCharacters("");
				w.flush();
				
				
				out.println(" var graph1_options = {");
				out.println("    numSequences: "+ variants.size() +",");
				out.println("    levels: [" + Arrays.stream(GenotypeType.values()).map(T->"\""+T.name()+"\"").collect(Collectors.joining(","))+"],");
				out.println("  }");

				out.println("  {");
				out.println("  var data=[");
				for(int i=0;i< samples.size();i++) {
					final String sn = samples.get(i);
					out.print("[");
					boolean first = true;
					for(int j=0;j < variants.size();j++ ) {
						final VariantContext ctx = variants.get(j);
						if(!ctx.isVariant()) continue;
						final Genotype gt = ctx.getGenotype(sn);
						if(gt.isHomRef()) continue;
						if(!first) out.print(",");
						out.println("{\"type\": \""+ ctx.getContig()+":"+ctx.getStart() +"\",\"level\": \""+ gt.getType().name()+"\"}");
						first=false;
						}
					out.print("]"+(i+1 < variants.size()?",":""));
					}
				out.println("    ];");
				out.println("console.log(data);");
				out.println("    graph_sandbox = new SequenceBraiding(data, '"+svgid+"', graph1_options)");

				
				out.println("    d3.selectAll('.lvlname').select('text').attr('x', 200)");
				out.println("    d3.selectAll('.lvlname').select('circle').attr('cx', 193)");
				out.println("    d3.selectAll('.path_top_text').text((d, i) => {");
				out.println("    switch(i) {");
				for(int i=0;i< samples.size();i++) {
					out.println("case "+i+": return \""+samples.get(i)+"\"; break;");
					}
				out.println("    }})");

				
				out.println("  }");
	
	
				
				
				out.flush();
				w.writeEndElement();//script
				w.writeEndElement();//html
	
				w.writeEndDocument();
				w.flush();
				out.flush();
				}
			return 0;
			}
		catch(Throwable err ) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new VcfBraiding().instanceMainWithExit(args);
	}
}
