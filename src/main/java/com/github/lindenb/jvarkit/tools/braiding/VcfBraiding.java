/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.braiding;

import java.io.PrintStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC
## Example

```
bcftools view src/test/resources/rotavirus_rf.vcf.gz "RF02" "RF03" |\
	java -jar /home/lindenb/src/jvarkit-git/dist/vcfbraiding.jar --title "Rotavirus Variants" > variants.html
```

## Screenshots

  https://twitter.com/yokofakun/status/1319221221611941889
  
  ![https://twitter.com/yokofakun/status/1319221221611941889](https://pbs.twimg.com/media/Ek7QRz9WMAApITu?format=jpg&name=large)

  https://twitter.com/yokofakun/status/1319228442043387905
  
  ![https://twitter.com/yokofakun/status/1319228442043387905](https://pbs.twimg.com/media/Ek7W9jaXEAEQ_TN?format=png&name=small)



## See also:
 
  * https://visdunneright.github.io/sequence_braiding/docs/

END_DOC

 */
@Program(name="vcfbraiding",
	description="visualization for variants and attributes using https://visdunneright.github.io/sequence_braiding/docs/ .",
	keywords={"vcf","visualization"},
	creationDate="20201021",
	modificationDate="20201022"
	)
public class VcfBraiding extends Launcher {
	private static final Logger LOG = Logger.build(VcfBraiding.class).make();
	private static final String DEFAULT_BASE = "https://visdunneright.github.io/sequence_braiding/docs/";
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-B","--base"},description="Base URL for code :" + DEFAULT_BASE)
	private String base = DEFAULT_BASE;
	@Parameter(names={"--id"},description="id for svg element.")
	private String svg_element_id = "vcfid";
	@DynamicParameter(names = "-D", description = "Dynamic parameters for API 'options'. -Dkey=value . Keys are currently: show_seq_names forceLevelName animate  colorbysequence width height fontSize padding.")
	private Map<String, String> parameters = new HashMap<>();
	@Parameter(names={"-T","--title"},description="title")
	private String title="";
	@Parameter(names={"--hom-ref","-hr"},description="remove sample that are all HOM_REF for all variants.")
	private boolean remove_all_hom_ref = false;
	@Parameter(names={"--no-call","-nc"},description="remove sample that are all NO_CALL for all variants.")
	private boolean remove_all_nocall = false;
	@Parameter(names={"--alleles"},description="show alleles in header.")
	private boolean show_alleles = false;



	private static String quote(String s) {
		return "\""+ s +"\"";
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final List<String> samples = new ArrayList<>();
			final List<VariantContext> variants = new ArrayList<>();
			final String input = oneFileOrNull(args);

			try(VCFIterator r= super.openVCFIterator(input)) {
				final VCFHeader header = r.getHeader();
				samples.addAll(header.getGenotypeSamples());
				while(r.hasNext()) {
					final VariantContext vc = r.next();
					if(!vc.isVariant()) continue;
					variants.add(vc);
					}
				}

			if(remove_all_hom_ref) samples.removeIf(S-> variants.stream().map(V->V.getGenotype(S)).allMatch(G->G.isHomRef()));
			if(remove_all_nocall) samples.removeIf(S-> variants.stream().map(V->V.getGenotype(S)).allMatch(G->G.isNoCall()));

			
			
			if(variants.size()<3) {
				LOG.error("Not enough variant was found. N=" + variants.size());
				return -1;
				}
			
			if(samples.size()<1) {
				LOG.error("Not enough samples was found. N=" + samples.size());
				return -1;
				}
			
			final boolean one_contig = variants.stream().map(V->V.getContig()).collect(Collectors.toSet()).size()==1;
			final Function<VariantContext, String> variantToStr = (V)->{
				final String s = String.valueOf(V.getStart()) + (show_alleles?":"+V.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining("/")):"");
				if(one_contig) return s;
				return V.getContig()+":"+s;
				};
			
			
			
			try(PrintStream out = super.openPathOrStdoutAsPrintStream(this.outputFile)) {
				
				final XMLOutputFactory factory= XMLOutputFactory.newInstance();
				final XMLStreamWriter w= factory.createXMLStreamWriter(out,"UTF-8");
				//w.writeStartDocument("UTF-8","1.0");
				w.writeStartElement("html");
				w.writeAttribute("lang","en");
				w.writeCharacters("\n");
				w.writeStartElement("head");
					w.writeStartElement("title");
						w.writeCharacters(StringUtil.isBlank(title)?(input==null?"Sequence Braiding":input):title);
					w.writeEndElement();//title
					
					w.writeEmptyElement("meta");
					w.writeAttribute("charset","utf-8");
	
					
					w.writeEmptyElement("meta");
					w.writeAttribute("name","description");
					w.writeAttribute("content",StringUtil.isBlank(title)?(input==null?"Sequence Braiding":input):title);
			
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
					w.writeCharacters("");
					w.writeEndElement();
			
					w.writeStartElement("link");
					w.writeAttribute("rel","stylesheet");
					w.writeAttribute("href",base + "css/skeleton.css");
					w.writeCharacters("");
					w.writeEndElement();
					
					w.writeStartElement("link");
					w.writeAttribute("rel","icon");
					w.writeAttribute("type","image/png");
					w.writeAttribute("href",base + "images/favicon.png");
					w.writeEndElement();
			
					w.writeStartElement("script");
					w.writeAttribute("src","https://d3js.org/d3.v5.min.js");
					w.writeCharacters("");
					w.writeEndElement();
					
				w.writeEndElement();//head
				w.writeStartElement("body");
				w.writeStartElement("h3");
				w.writeCharacters(StringUtil.isBlank(this.title)?(getClass().getSimpleName()+(input==null?"":" : " +input)):title);
				w.writeEndElement();
				
	
				w.writeStartElement("div");
				w.writeStartElement("svg");
				w.writeAttribute("id",svg_element_id);
				w.writeAttribute("style","margin-left: 5%");
				w.writeEndElement();//svg
				w.writeEndElement();//div
				w.writeEmptyElement("hr");
				
				w.writeStartElement("div");
				w.writeCharacters("Made with "+getClass().getSimpleName()+" version:"+getVersion()+". Pierre Lindenbaum PhD. 2020. See also:");
				w.writeStartElement("a");
				w.writeAttribute("href", "https://visdunneright.github.io/sequence_braiding/docs/");
				w.writeCharacters("https://visdunneright.github.io/sequence_braiding/docs/");
				w.writeEndElement();//a
				w.writeCharacters(".");
				w.writeEndElement();//div
				w.writeEndElement();//body
				
				
				for(String js : new String[] {
						"../lib/align_pair_quad.js",
						"../lib/pairwise_alignment_dna.js",
						"../js/util.js",
						"../js/SequenceBraiding.js",
						}) {
					w.writeStartElement("script");
					w.writeAttribute("src",base+js);
					w.writeCharacters("");
					w.writeEndElement();//script
					}
				//
				w.writeStartElement("script");
				w.writeAttribute("type","text/javascript");
				w.writeCharacters("\n");
				w.flush();
				
				final boolean show_seq_names = this.parameters.getOrDefault("show_seq_names", "true").equals("true");
				
				out.println("var options = {");
				out.println("   numSequences: "+samples.size()+",");
				out.println("   show_seq_names: "+show_seq_names+",");
				out.println("   forceLevelName: "+this.parameters.getOrDefault("forceLevelName", "true")+",");
				out.println("   animate: "+this.parameters.getOrDefault("animate", "false")+",");
				out.println("   colorbysequence: "+this.parameters.getOrDefault("colorbysequence", "true")+",");
				if(this.parameters.containsKey("fontSize")) out.println("   fontSize: "+quote(this.parameters.getOrDefault("fontSize", "0.9em"))+",");
				out.println("   width: "+this.parameters.getOrDefault("width", "window.innerWidth*0.9")+",");
				out.println("   height: "+this.parameters.getOrDefault("height", "800")+",");
				if(this.parameters.containsKey("padding"))  out.println("   padding: "+ this.parameters.getOrDefault("padding", "200")+",");
				out.println("   levels: [ "+ quote("NO_CALL") +"," + quote("HOM_REF")+","+ quote("HET")+","+quote("HOM_VAR")+"]");
				out.println("  }");

				out.println("  {");
				out.println("  var data=[");
				for(int i=0;i< samples.size();i++) {
					final String sn = samples.get(i);
					out.print((i==0?"":",") + "[");
					out.print(variants.stream().
						map(V->"\t{\"type\": "+ quote(variantToStr.apply(V)) +",\"level\": "+ quote(V.getGenotype(sn).getType().name())+
								(show_seq_names?","+quote("seq_name")+":"+quote(sn):"")+"}").
						collect(Collectors.joining(",")));
					out.print("]\n");
					}
				out.println("];");
				//out.println("console.log(data);");
				
				out.println("    graph_sandbox = new SequenceBraiding(data, '"+svg_element_id+"', options);");

				
				//out.println("    d3.selectAll('.lvlname').select('text').attr('x', 200)");
				//out.println("    d3.selectAll('.lvlname').select('circle').attr('cx', 193)");
				out.println("    d3.selectAll('.path_top_text').text((d, i) => {");
				out.println("    switch(i) {");
				for(int i=0;i< variants.size();i++) {
					out.println("case "+i+": return "+quote(variantToStr.apply(variants.get(i)))+"; break;");
					}
				out.println("    }});");

				
				out.println("  }");
	
				
				out.flush();
				w.writeEndElement();//script
				w.writeEndElement();//html
	
				//w.writeEndDocument();
				w.flush();
				out.flush();
				}
			return 0;
			}
		catch(final Throwable err ) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new VcfBraiding().instanceMainWithExit(args);
	}
}
