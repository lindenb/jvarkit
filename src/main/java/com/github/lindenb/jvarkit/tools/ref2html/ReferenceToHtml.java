/*

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
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.ref2html;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.lang.StaticCodeExtractor;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;
import com.google.gson.stream.JsonWriter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/*
BEGIN_DOC

## Example

```
$ java -jar dist/ref2html.jar -o TMP --regions "RF02:1-500" -R src/test/resources/rotavirus_rf.fa src/test/resources/rotavirus_rf.vcf.gz
```

END_DOC
*/

@Program(name="ref2html",
description="Produce a html view of eah interval from a VCF file",
keywords={"vcf","html"},
creationDate="20211012",
modificationDate="20211013",
generate_doc=false
)
public class ReferenceToHtml extends Launcher {
	private static final Logger LOG = Logger.build(ReferenceToHtml.class).make();

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path archiveOutput = null ;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null ;
	@Parameter(names={"-r","--region","--regions","--interval"},description=IntervalListProvider.OPT_DESC,required=true)
	private String regionsStr = null;
	@Parameter(names={"---prefix"},description="file prefix.")
	private String prefix = "";

	
	private void writeMeta(final XMLStreamWriter w) throws XMLStreamException {
	w.writeEmptyElement("meta");
		w.writeAttribute("http-equiv", "Content-Type");
		w.writeAttribute("content", "text/html; charset=utf-8");
	w.writeEmptyElement("meta");
		w.writeAttribute("http-equiv", "author");
		w.writeAttribute("content", "Pierre Lindenbaum");
	}
	
	private void checkBox(final XMLStreamWriter out,String id,String label) throws XMLStreamException{
	out.writeEmptyElement("input");
	out.writeAttribute("id", id);
	out.writeAttribute("type", "checkbox");
	out.writeAttribute("checked", "true");
	out.writeStartElement("label");
	out.writeAttribute("for",id);
	out.writeCharacters(label);
	out.writeEndElement();
	}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			try(VCFReader reader  = VCFReaderFactory.makeDefault().open(Paths.get(oneAndOnlyOneFile(args)),true)) {
				final VCFHeader header= reader.getHeader();
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
				final Optional<String> buildName = SequenceDictionaryUtils.getBuildName(dict);
				try( ReferenceSequenceFile reference = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.faidx)) {
					SequenceUtil.assertSequenceDictionariesEqual(SequenceDictionaryUtils.extractRequired(reference),dict);
					final List<Locatable> locs = IntervalListProvider.from(this.regionsStr).dictionary(dict).stream().
						sorted(new ContigDictComparator(dict).createLocatableComparator()).
						collect(Collectors.toList());
					if(locs.isEmpty()) {
						LOG.error("no region was defined");
						return -1;
						}
					
					final XMLOutputFactory xof = XMLOutputFactory.newFactory();
					
					try(ArchiveFactory archive = ArchiveFactory.open(archiveOutput)) {
						try(PrintWriter pw = archive.openWriter(this.prefix + "script.js")) {
							pw.println( StaticCodeExtractor.forClass(this.getClass()).extract("SCRIPT").get());
							pw.flush();
						}
						try(PrintWriter pw = archive.openWriter(this.prefix + "style.css")) {
							pw.println( StaticCodeExtractor.forClass(this.getClass()).extract("CSS").get());
	
							pw.flush();
							}
					
						final OutputStream index_os = archive.openOuputStream(this.prefix + "index.html");
						final XMLStreamWriter index_html = xof.createXMLStreamWriter(index_os, "UTF-8");
						index_html.writeStartDocument("UTF-8", "1.0");
						index_html.writeStartElement("html");
						
						index_html.writeStartElement("head");
						writeMeta(index_html);											
						index_html.writeStartElement("title");
						index_html.writeCharacters(this.getProgramName());
						index_html.writeEndElement();//title
						
						index_html.writeStartElement("link");
						index_html.writeAttribute("rel", "stylesheet");
						index_html.writeAttribute("href", this.prefix +  "style.css");
						index_html.writeCharacters("");
						index_html.writeEndElement();//link

						index_html.writeEndElement();//head
						
						index_html.writeStartElement("body");
						index_html.writeStartElement("ul");
						
						
						for(final Locatable loc:locs) {
							final List<VariantContext> variants;
							try(CloseableIterator<VariantContext> iter = reader.query(loc)) {
								variants = iter.stream().collect(Collectors.toList());
								}
							
							StringWriter sw = new StringWriter();
							
							sw.append("var g = ");
							JsonWriter jsw = new JsonWriter(sw);
							jsw.beginObject();
							jsw.name("contig").value(loc.getContig());
							jsw.name("start").value(loc.getStart());
							jsw.name("end").value(loc.getEnd());
							jsw.name("length").value(loc.getLengthOnReference());
							jsw.name("sequence").value(reference.getSubsequenceAt(loc.getContig(), loc.getStart(), loc.getEnd()).getBaseString());
							jsw.name("variants");
							jsw.beginArray();
							for(VariantContext ctx:variants) {
								jsw.beginObject();
								jsw.name("start").value(ctx.getStart());
								jsw.name("end").value(ctx.getEnd());
								jsw.name("id").value(ctx.hasID()?ctx.getID():"");
								jsw.name("type").value(ctx.getType().name());
								jsw.name("ref").value(ctx.getReference().getDisplayString());
								jsw.name("alts");
								jsw.beginArray();
								for(Allele alt:ctx.getAlternateAlleles()) {
									jsw.value(alt.getDisplayString());
									}
								jsw.endArray();
								jsw.endObject();
								}
							jsw.endArray();
							jsw.endObject();
							jsw.flush();
							jsw.close();
							sw.append(";");
							
							
							
							
							
							final String filename = this.prefix + loc.getContig()+"_"+loc.getStart()+"_"+loc.getEnd()+".html";

							final String title= (buildName.isPresent()?buildName.get()+" ":"")+new SimpleInterval(loc).toNiceString();
							OutputStream os = archive.openOuputStream(filename);
							XMLStreamWriter out = xof.createXMLStreamWriter(os, "UTF-8");
							out.writeStartDocument("UTF-8", "1.0");
							out.writeStartElement("html");
							out.writeStartElement("head");
							writeMeta(out);
							out.writeStartElement("script");
							out.writeCharacters(sw.toString());	
							out.writeEndElement();//script
							out.writeStartElement("script");
							out.writeAttribute("src", this.prefix +  "script.js");
							out.writeCharacters("");
							out.writeEndElement();//script
							
							out.writeStartElement("link");
							out.writeAttribute("rel", "stylesheet");
							out.writeAttribute("href", this.prefix +  "style.css");
							out.writeCharacters("");
							out.writeEndElement();//link

							
							out.writeStartElement("title");
							out.writeCharacters(title+" N="+StringUtils.niceInt(variants.size()));
							out.writeEndElement();//title
							out.writeEndElement();//head
							out.writeStartElement("body");
							
							out.writeStartElement("h1");
							out.writeCharacters(title+" N="+StringUtils.niceInt(variants.size()));
							out.writeEndElement();
							
							out.writeStartElement("div");
							checkBox(out,"showPos","Line Number");
							checkBox(out,"showSpace","Spaces");
							checkBox(out,"showVar","Show Variants");
							out.writeCharacters(" ");
							out.writeEmptyElement("input");
								out.writeAttribute("id", "primer");
								out.writeAttribute("type", "text");
								out.writeAttribute("placeholder", "Primer Sequence");
							out.writeStartElement("button");
								out.writeAttribute("id", "search");
								out.writeCharacters("Search...");
							out.writeEndElement();
							out.writeEndElement();//div
							
							
							out.writeStartElement("div");
							out.writeAttribute("class","sequence");
							out.writeAttribute("id","main");
							out.writeEndElement();//div
							
							out.writeEndElement();//body
							out.writeEndElement();//html
							os.flush();
							os.close();
							
							index_html.writeStartElement("li");
							index_html.writeStartElement("a");
							index_html.writeAttribute("href", filename);
							index_html.writeCharacters(title+" N="+StringUtils.niceInt(variants.size()));
							index_html.writeEndElement();//a
							index_html.writeEndElement();//li
							}
						index_html.writeEndElement();//ul
						index_html.writeEndElement();//body
						index_html.writeEndElement();//html
						index_html.writeEndDocument();
						index_html.close();
						index_os.flush();
						index_os.close();
						}
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new ReferenceToHtml().instanceMainWithExit(args);

	}

}

/**

BEGIN_SCRIPT

var g2={
	"fastalen":60,
	"hits":[]
	};
function pad(num,len) {
	var s=""+num;
	while(s.length<len) s=" "+s;
	return s;
	}
function removeAllChildNodes(parent) {
    while (parent.firstChild)  parent.removeChild(parent.firstChild);
    return parent;
	}
function revcomp(s) {
	var c,s2="",i;
	for(i=0;i<s.length;i++) {
		switch(s[i]) {
			case ' ' : break;
			case 'a':case 'A': c = 'T'; break;
			case 't':case 'T': c = 'A'; break;
			case 'g':case 'G': c = 'C'; break;
			case 'c':case 'C': c = 'G'; break;
			default: c= 'N'; break;
			}
		s2 = c+s2;
		}
	return s2;
	}
function findHits(query) {
	var hits=[];
	var i,j,side,L1 = query.length;
	query=query.toUpperCase();
	if(query.length>0) {
		for(side=0;side<2;side++) {
			var s;
			if(side==0) {
				s = query;
				} else {
				s= revcomp(query);
				if(s==query) break;
				}
			
			for(i=0;i+L1<=g.length;i++) {
				for(j=0;j< L1;j++) {
					if(s[j]!=g.sequence[i+j].toUpperCase()) break;
					}
				if(j==L1) {
					hits.push({"start":g.start+i,"end": g.start+(i+L1-1),"strand":side});
					}
				}
			}
		}
	return hits;
	}

function findPrimer()  {
	var dna = document.getElementById("primer").value.trim().toUpperCase();
	if(!dna.match(/^[ATGC]+$/))
		{
		g2.hits= [];
		} else {
		g2.hits=findHits(dna);
		}
	drawSequence();
	}

function findVariantsAt(pos) {
	var i,array=[];
	for(i=0;i< g.variants.length;i++) {
		var ctx = g.variants[i];
		if(ctx.start == pos) array.push(ctx);
		}
	return array;
	}
function findHitsAt(pos) {
	var i,array=[];
	for(i=0;i< g2.hits.length;i++) {
		var p = g2.hits[i];
		if(p.start<= pos && pos<= p.end) array.push(p);
		}
	return array;
	}

function variantToString(ctx) {
	return ""+ ctx.start+ " "+ctx.type +" "+ctx.ref+"/"+ctx.alts.join("/")+(ctx.id==""?"":" ("+ctx.id+")");
	}
	
function drawSequence() {
	var i,j;
	var pre = removeAllChildNodes(document.getElementById("main"));
	var linenumber = document.getElementById("showPos").checked;
	var showspace = document.getElementById("showSpace").checked;
	var showVar = document.getElementById("showVar").checked;
	console.log(""+linenumber);
	for(i=0;i< g.length;i++) {
		if(i%g2.fastalen==0) {
			if(i>0) {
				pre.appendChild(document.createTextNode("\n"));

				}
			if(linenumber) pre.appendChild(document.createTextNode(pad(g.start + i,9)+"  "));
		} else if(i%10==0 && showspace) {
			pre.appendChild(document.createTextNode(" "));
			}
		
		var variants = findVariantsAt(g.start+i);
		var hits = findHitsAt(g.start+i);
		var clazz=(variants.length>0?"variant":"")+" "+(hits.length>0?"primer":"");
		clazz=clazz.trim();
		
		if(clazz.trim()!="") {
			var span = document.createElement("span");
			span.setAttribute("class",clazz);
			if(showVar && variants.length>0) {
				for(j=0;j < variants.length;j++) {
					span.appendChild(document.createTextNode("["+variantToString(variants[j])+"]"));
					}
				}
			else {
				span.appendChild(document.createTextNode(g.sequence[i]));	
				}
			pre.appendChild(span);
			}
		else
			{
			pre.appendChild(document.createTextNode(g.sequence[i]));
			}
		}
	}

function docLoaded() {
	drawSequence();
	var t  = document.getElementById("showPos");
	t.addEventListener("change",drawSequence,false);
	t  = document.getElementById("showSpace");
	t.addEventListener("change",drawSequence,false);
	t  = document.getElementById("showVar");
	t.addEventListener("change",drawSequence,false);
	t = document.getElementById("search");
	t.addEventListener("click",findPrimer,false);
}

window.addEventListener("load",docLoaded,false)


END_SCRIPT


BEGIN_CSS
.variant {
	color:red;
	}
.primer {
	text-decoration: underline;
	}

div.sequence {   
	display: block;
    unicode-bidi: embed;
    font-family: monospace;
    white-space: pre;
    background-color:lightgray;
    }
END_CSS


**/
