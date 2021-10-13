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
package com.github.lindenb.jvarkit.tools.vcf2html;

import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
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
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;

@Program(name="vcf2html",
description="Produce a html view of eah interval from a VCF file",
keywords={"vcf","html"},
creationDate="20211012",
modificationDate="20211012",
generate_doc=false
)
public class VcfToHtml extends Launcher {
	private static final Logger LOG = Logger.build(VcfToHtml.class).make();

	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC,required=true)
	private Path archiveOutput = null ;

	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx = null ;
	@Parameter(names={"-r","--region","--regions","--interval"},description=IntervalListProvider.OPT_DESC,required=true)
	private String regionsStr = null;
	
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
						final OutputStream index_os = archive.openOuputStream("index.html");
						XMLStreamWriter index_html = xof.createXMLStreamWriter(index_os, "UTF-8");
						index_html.writeStartDocument("UTF-8", "1.0");
						index_html.writeStartElement("html");
						index_html.writeStartElement("body");
						
						
						
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
							
							
							
							
							
							final String filename = loc.getContig()+"_"+loc.getStart()+"_"+loc.getEnd()+".html";

							final String title= (buildName.isPresent()?buildName.get()+" ":"")+new SimpleInterval(loc).toNiceString();
							OutputStream os = archive.openOuputStream(filename);
							XMLStreamWriter out = xof.createXMLStreamWriter(os, "UTF-8");
							out.writeStartDocument("UTF-8", "1.0");
							out.writeStartElement("html");
							out.writeStartElement("head");
							out.writeStartElement("script");
							out.writeCharacters(sw.toString());	
							
							out.writeCharacters(
									"function subref(start,end) {}" +
									"function createCtx(ith) {}" +
									"function repaint() {}"
								);
							
							out.writeEndElement();//script
							out.writeStartElement("script");
							out.writeAttribute("src", "script.js");
							out.writeEndElement();//script
							out.writeStartElement("title");
							out.writeCharacters(title+" N="+StringUtils.niceInt(variants.size()));
							out.writeEndElement();//title
							out.writeEndElement();//head
							out.writeStartElement("body");
							
							out.writeEmptyElement("input");
							
							
							out.writeStartElement("div");
							out.writeEndElement();//div
							
							out.writeEndElement();//body
							out.writeEndElement();//html
							os.flush();
							os.close();
							
							index_html.writeStartElement("li");
							index_html.writeStartElement("a");
							index_html.writeAttribute("a", filename);
							index_html.writeCharacters(title+" N="+StringUtils.niceInt(variants.size()));
							index_html.writeEndElement();//a
							index_html.writeEndElement();//li
							}
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
	
	public static void main(String[] args) {
		new VcfToHtml().instanceMainWithExit(args);

	}

}
