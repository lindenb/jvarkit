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
package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.LineIterator;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/*
BEGIN_DOC


## Example

```bash
$ java -jar dist/gtf2xml.jar src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | xmllint --format - | head -n 100
<?xml version="1.0" encoding="UTF-8"?>
<gtf genebuild-last-updated="2013-09" genome-build="GRCh37.p13" genome-build-accession="NCBI:GCA_000001405.14" genome-date="2009-02" genome-version="GRCh37">
  <gene id="ENSG00000100403" chrom="22" start="41697526" end="41756151" strand="+" source="ensembl_havana" type="gene">
    <attributes>
      <attribute key="gene_id">ENSG00000100403</attribute>
      <attribute key="gene_version">10</attribute>
      <attribute key="gene_name">ZC3H7B</attribute>
      <attribute key="gene_source">ensembl_havana</attribute>
      <attribute key="gene_biotype">protein_coding</attribute>
    </attributes>
    <transcripts>
      <transcript id="ENST00000486331" chrom="22" start="41697719" end="41732847" strand="+" source="havana" type="transcript">
        <attributes>
          <attribute key="gene_id">ENSG00000100403</attribute>
          <attribute key="gene_version">10</attribute>
          <attribute key="transcript_id">ENST00000486331</attribute>
          <attribute key="transcript_version">1</attribute>
          <attribute key="gene_name">ZC3H7B</attribute>
          <attribute key="gene_source">ensembl_havana</attribute>
          <attribute key="gene_biotype">protein_coding</attribute>
          <attribute key="transcript_name">ZC3H7B-002</attribute>
          <attribute key="transcript_source">havana</attribute>
          <attribute key="transcript_biotype">retained_intron</attribute>
          <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
          <attribute key="havana_transcript_version">1</attribute>
        </attributes>
        <exon chrom="22" start="41697719" end="41697776" strand="+" source="havana" type="exon">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="exon_number">1</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
            <attribute key="exon_id">ENSE00001942555</attribute>
            <attribute key="exon_version">1</attribute>
          </attributes>
        </exon>
        <transcript chrom="22" start="41697719" end="41732847" strand="+" source="havana" type="transcript">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
          </attributes>
        </transcript>
        <exon chrom="22" start="41716659" end="41716717" strand="+" source="havana" type="exon">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="exon_number">2</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
            <attribute key="exon_id">ENSE00003530265</attribute>
            <attribute key="exon_version">1</attribute>
          </attributes>
        </exon>
        <exon chrom="22" start="41721568" end="41721601" strand="+" source="havana" type="exon">
          <attributes>
            <attribute key="gene_id">ENSG00000100403</attribute>
            <attribute key="gene_version">10</attribute>
            <attribute key="transcript_id">ENST00000486331</attribute>
            <attribute key="transcript_version">1</attribute>
            <attribute key="exon_number">3</attribute>
            <attribute key="gene_name">ZC3H7B</attribute>
            <attribute key="gene_source">ensembl_havana</attribute>
            <attribute key="gene_biotype">protein_coding</attribute>
            <attribute key="transcript_name">ZC3H7B-002</attribute>
            <attribute key="transcript_source">havana</attribute>
            <attribute key="transcript_biotype">retained_intron</attribute>
            <attribute key="havana_transcript">OTTHUMT00000320697</attribute>
            <attribute key="havana_transcript_version">1</attribute>
            <attribute key="exon_id">ENSE00003553644</attribute>
            <attribute key="exon_version">1</attribute>
          </attributes>
        </exon>
        (...)
```

END_DOC
 */
@Program(name="gtf2xml",
	description="Convert GTF/GFF to XML",
	keywords={"xml","gtf","gff","gff3"},
	creationDate="20150811",
	modificationDate="20190823",
	biostars={478242}
	)
public class Gtf2Xml extends Launcher{
	private static final Logger LOG = Logger.build(FixVCF.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-f","--features"},description="Don't record features types.")
	private boolean disable_feature_type=false;
	@Parameter(names={"-s","--sources"},description="Don't record sources")
	private boolean disable_sources=false;
	@Parameter(names={"-d","--dict"},description="Don't record contigs")
	private boolean disable_dict=false;
	@Parameter(names={"-a","--attributes"},description="Don't record attribute types.")
	private boolean disable_att_keys=false;
	
	private final Map<String,Long> seqdict=new LinkedHashMap<>();
	private final Set<String> att_keys=new HashSet<>();
	private final Set<String> sources=new HashSet<>();
	private final Set<String> types=new HashSet<>();
	
	
	private void commonFields(final XMLStreamWriter w,final GTFLine line) throws XMLStreamException
		{
		w.writeAttribute("chrom", line.getContig());
		w.writeAttribute("start",String.valueOf(line.getStart()));
		w.writeAttribute("end",String.valueOf(line.getEnd()));

		if(!this.disable_dict) {
			Long contifLength  = this.seqdict.get(line.getContig());
			if(contifLength==null) contifLength=0L;
			this.seqdict.put(line.getContig(),Math.max(line.getEnd(), contifLength));
			}
		if(line.getScore()!=null)  {
			w.writeAttribute("score",String.valueOf(line.getScore()));
		}
		if(line.getStrand()!=GTFLine.NO_STRAND) {
			w.writeAttribute("strand",String.valueOf(line.getStrand()));
		}
		if(!StringUtil.isBlank(line.getSource()))  {
			if(!this.disable_sources) this.sources.add(line.getSource());
			w.writeAttribute("source",line.getSource());
		}
		if(!StringUtil.isBlank(line.getType()))  {
			if(!disable_feature_type) this.types.add(line.getType());
			w.writeAttribute("type",line.getType());
		}
		if(line.getPhase()!=GTFLine.NO_PHASE) {
			w.writeAttribute("phase",String.valueOf(line.getPhase()));
			}
		
		w.writeStartElement("attributes");
		for(final Iterator<Map.Entry<String, String>> kvr=line.getAttributeIterator();
				kvr.hasNext();
				)
			{
			final Map.Entry<String, String> kv = kvr.next();
			w.writeStartElement("attribute");
			w.writeAttribute("key", kv.getKey());
			if(!this.disable_att_keys) this.att_keys.add(kv.getKey());
			w.writeCharacters(kv.getValue());
			w.writeEndElement();
			}
		w.writeEndElement();
			
		}
	
	private void writeTranscript(final XMLStreamWriter w,final String transcript_id,final List<GTFLine> lines) throws XMLStreamException,IOException
	{
	final Optional<GTFLine> optGeneLine = lines.stream().filter(G->G.getType().equals("transcript")).findFirst();
	if(!optGeneLine.isPresent()) throw new RuntimeIOException("No 'transcript' found for transcript : "+transcript_id+" in "+lines);
	w.writeStartElement("transcript");
	w.writeAttribute("id", transcript_id);
	commonFields(w,optGeneLine.get());
	for(final GTFLine l:lines) {
		write(w, l);
		}
	
	w.writeEndElement();
	}
	
	private void writeGene(final XMLStreamWriter w,final String gene_id,final List<GTFLine> lines) throws XMLStreamException,IOException
		{
		final Optional<GTFLine> optGeneLine = lines.stream().filter(G->G.getType().equals("gene")).findFirst();
		if(!optGeneLine.isPresent()) throw new RuntimeIOException("No 'gene' found for gene_id : "+gene_id);
		w.writeStartElement("gene");
		w.writeAttribute("id", gene_id);
		commonFields(w,optGeneLine.get());
		
		final Map<String,List<GTFLine>> id2transcript = new HashMap<>();
		for(final GTFLine line:lines) {
			final String transcript_id = line.getAttribute("transcript_id");
			if(StringUtil.isBlank(transcript_id))  continue;
			
			List<GTFLine> lines2 = id2transcript.get(transcript_id);
			if(lines2==null) {
				lines2 = new ArrayList<>();
				id2transcript.put(transcript_id,lines2);
				}
			lines2.add(line);
			}
		if(!id2transcript.isEmpty()) {
			w.writeStartElement("transcripts");
			for(final String transcript_id:id2transcript.keySet()) {
				writeTranscript(w,transcript_id,id2transcript.get(transcript_id));
				}
			w.writeEndElement();
			}
		
		w.writeEndElement();
		}

	
	private void write(final XMLStreamWriter w,final GTFLine line) throws XMLStreamException,IOException
		{
		w.writeStartElement(line.getType());
		
		commonFields(w,line);
		
		
		w.writeEndElement();
		w.writeCharacters("\n");
		}
	
	@Override
	public int doWork(final List<String> args) {

		LineIterator r=null;
		XMLStreamWriter w=null;
		PrintWriter fw=null;
		try {
			String inputName=oneFileOrNull(args);
			r = (StringUtil.isBlank(inputName)?
					IOUtils.openStreamForLineIterator(stdin()):
					IOUtils.openURIForLineIterator(inputName)
					);
			XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w = xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				w = xof.createXMLStreamWriter((fw=super.openPathOrStdoutAsPrintWriter(outputFile)));
				}
			final GTFCodec codec = new GTFCodec();
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("gtf");
			
			final GTFCodec.GTFHeader header = codec.readActualHeader(r);
			for(final String headerLine:header.getLines()) {
				if(!headerLine.startsWith("#!")) continue;
				final int ws = headerLine.indexOf(' ');
				if(ws==-1) continue; //??
				w.writeAttribute(
						headerLine.substring(2, ws),
						headerLine.substring(ws+1).trim());
				
				}
		
			final Map<String,List<GTFLine>> geneid2lines = new HashMap<>(50_000);
			
			while(r.hasNext())
				{
				final String line=r.next();
				final GTFLine gtfline = codec.decode(line);
				if(gtfline==null) continue;
				final String gene_id = gtfline.getAttribute("gene_id");
				if(StringUtil.isBlank(gene_id)) {
					write(w,gtfline);
					}
				else
					{
					List<GTFLine> lines = geneid2lines.get(gene_id);
					if(lines==null) {
						lines = new ArrayList<>();
						geneid2lines.put(gene_id,lines);
						}
					lines.add(gtfline);
					}
				}
			
			for(final String gene_id:geneid2lines.keySet()) {
				writeGene(w,gene_id,geneid2lines.get(gene_id));
				}
			
			if(!this.disable_att_keys) {
				w.writeStartElement("attributes");
				for(final String k : this.att_keys)
					{
					w.writeStartElement("attribute");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			
			
			if(!this.disable_feature_type) {
				w.writeStartElement("types");
				for(String k : this.types)
					{
					w.writeStartElement("type");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}
			
			if(!this.disable_sources) {
				w.writeStartElement("sources");
				for(final String k : this.sources)
					{
					w.writeStartElement("source");
					w.writeCharacters(k);
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}

			
			if(!this.disable_dict) {
				w.writeStartElement("dict");
				for(final String k : this.seqdict.keySet())
					{
					w.writeEmptyElement("chrom");
					w.writeAttribute("name",k);
					w.writeAttribute("length", String.valueOf(this.seqdict.get(k)));
					w.writeCharacters("\n");
					}
				w.writeEndElement();
				w.writeCharacters("\n");
				}

			
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			return 0;
			}
		catch (Exception e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(fw);
			}
		}

	
	public static void main(String[] args) throws IOException
		{
		new Gtf2Xml().instanceMainWithExit(args);
		}
}
