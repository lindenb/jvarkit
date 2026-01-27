/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gtf2xml;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.dict.DictionaryXmlSerializer;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.gtf.GTFCodec;
import com.github.lindenb.jvarkit.gtf.GTFLine;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.math.MinMaxInteger;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.samtools.util.Pileup;
import com.github.lindenb.jvarkit.tools.misc.FixVCF;
import com.github.lindenb.jvarkit.util.AutoMap;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
/*
BEGIN_DOC


## Example

(the xml in the example below is old)

```bash
$ java -jar dist/jvarkit.jar gtf2xml src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | xmllint --format - | head -n 100
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
	modificationDate="20260118",
	biostars={478242},
	jvarkit_amalgamion = true,
	menu="GTF/GFF Manipulation"
	)
public class Gtf2Xml extends Launcher{
	private static final Logger LOG = Logger.of(FixVCF.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--reference","-R"},description=DICTIONARY_SOURCE)
	private Path faidxPath = null;
	@Parameter(names={"--skip-attributes"},description="Don't print the following attributes. Multiple separated by commas/spaces/semicolons")
	private String skipAttributeStr = "";
	@Parameter(names={"--interval","--regions","--region","-r"},description="only in the following interval CHR:START-END. First Scan for genes in the interval and extend the interval to get the whole genes. Requires GTF input is a tabix index file.")
	private String intervalStr = "";
	@Parameter(names={"-d","--distance"},description="if >=0, add a 'y' attribute that could be used to display the bed records in a browser, " 
			+ "	this 'y' is the graphical row where the item should be displayed. " 
			+ "	This distance is the distance between two item where there is a collision. Memory consuming . "+
			DistanceParser.OPT_DESCRIPTION,
			converter = DistanceParser.StringConverter.class, splitter = NoSplitter.class)
	private int distance=-1;
	@Parameter(names={"--simple"},description="Don't print group data by gene/transcript. Print each GTF record on the fly")
	private boolean disable_group_by_gene = false;

	
	private final Set<String> skip_attributes_set = new HashSet<>();
	

	private static class CountRows {
		final int y;
		final int nrows;
		CountRows(int y, int nrows) {
			this.y = y;
			this.nrows = nrows;
			}
		}
	
	private void writeStartRecord(final XMLStreamWriter w,final GTFLine line,
			final CountRows countRows,
			final Locatable interval) throws XMLStreamException
		{
		w.writeStartElement(line.getType().toLowerCase());
		w.writeAttribute("chrom", line.getContig());
		w.writeAttribute("start",String.valueOf(line.getStart()));
		w.writeAttribute("end",String.valueOf(line.getEnd()));
		if(interval!=null) {
			if(interval.overlaps(line)) {
				w.writeAttribute("f0",String.valueOf( (Math.max(interval.getStart(), line.getStart())-interval.getStart())/(double)interval.getLengthOnReference()));
				w.writeAttribute("f1",String.valueOf( (Math.min(interval.getEnd(), line.getEnd())-interval.getStart()+1)/(double)interval.getLengthOnReference()));
				}
			else
				{
				w.writeAttribute("visible","false");
				}
			}
		
		
		if(line.getScore()!=null)  {
			w.writeAttribute("score",String.valueOf(line.getScore()));
		}
		if(line.getStrand()!=GTFLine.NO_STRAND) {
			w.writeAttribute("strand",String.valueOf(line.getStrand()));
		}
		
		if(line.getPhase()!=GTFLine.NO_PHASE) {
			w.writeAttribute("phase",String.valueOf(line.getPhase()));
			}
		
		if(countRows!=null ) {
			w.writeAttribute("y",String.valueOf(countRows.y));
			w.writeAttribute("nrows",String.valueOf(countRows.nrows));
			}
	
		
		w.writeStartElement("attributes");
		final Map<String,String> atts = line.getAttributes();
		for(String key: atts.keySet()) {
			if(this.skip_attributes_set.contains(key)) continue;
			w.writeStartElement("attribute");
			w.writeAttribute("key", key);
			w.writeCharacters(atts.get(key));
			w.writeEndElement();
			}
		w.writeEndElement();//attributes
		}
	
	private void writeEndRecord(final XMLStreamWriter w,final GTFLine line) throws XMLStreamException
		{
		w.writeEndElement();
		w.writeCharacters("\n");
		}
	
	
	private void write(final XMLStreamWriter w,final GTFLine line, final CountRows countRows,final Locatable loc) throws XMLStreamException
		{
		writeStartRecord(w, line, countRows, loc);
		writeEndRecord(w, line);
		}
	
	
	
	private void compile(final XMLStreamWriter w, List<GTFLine> records,final Locatable interval) throws XMLStreamException {
		final List<GTFLine> genes = new ArrayList<>();
		final AutoMap<String, GTFLine, List<GTFLine>> gene2transcripts = AutoMap.makeList();
		final AutoMap<String, GTFLine, List<GTFLine>> transcript_components = AutoMap.makeList();
		
		while(!records.isEmpty()) {
			final GTFLine rec = records.remove(records.size()-1);
			if(rec.isGene()) {
				genes.add(rec);
				}
			else if(rec.isTranscript()) {
				final String gene_id = rec.getGeneId();
				if(StringUtils.isBlank(gene_id)) {
					LOG.warn("no gene_id for " + rec);
					continue;
					}
				gene2transcripts.insert(gene_id, rec);
				}
			else {
				final String transcript_id = rec.getTranscriptId();
				if(StringUtils.isBlank(transcript_id)) {
					LOG.warn("no transcript_id for " + rec);
					continue;
					}
				transcript_components.insert(transcript_id, rec);
				}
			}
		final Map<String,CountRows> gene2y = new HashMap<>();
		final Map<String,CountRows> transcript2y = new HashMap<>();
		if(this.distance>=0) {
			/* transcripts */
			final Pileup<GTFLine> pileup = new Pileup<>((LEFT,RIGHT)->CoordMath.getLength(LEFT.getEnd(),RIGHT.getStart())>= Gtf2Xml.this.distance);
			pileup.addAll(
				gene2transcripts.values()
					.stream()
					.flatMap(L->L.stream())
					.sorted((A,B)->Integer.compare(A.getStart(),B.getStart()))
					.collect(Collectors.toList())
					);
			int y=0;
			for(List<GTFLine> row : pileup.getRows()) {
				for(GTFLine rec:row) {
					final String transcript_id = rec.getTranscriptId();
					if(StringUtils.isBlank(transcript_id)) {
						throw new IllegalStateException();
						}
					transcript2y.put(transcript_id, new CountRows(y,pileup.getRowCount()));
					}
				++y;
				}
			pileup.clear();
			
			/* genes */
			Collections.sort(genes,(A,B)->Integer.compare(A.getStart(),B.getStart()));
			pileup.addAll(genes);
			y=0;
			for(List<GTFLine> row : pileup.getRows()) {
				for(GTFLine rec:row) {
					final String gene_id = rec.getGeneId();
					if(StringUtils.isBlank(gene_id)) {
						throw new IllegalStateException();
						}
					gene2y.put(gene_id,new CountRows(y,pileup.getRowCount()));
					}
				++y;
				}
			}
		
		for(GTFLine gene: genes) {
			final String gene_id = gene.getGeneId();
			writeStartRecord(w, gene, gene2y.get(gene_id), interval);

			final List<GTFLine> transcripts = gene2transcripts.getOrDefault(gene_id,Collections.emptyList());
			w.writeStartElement("transcripts");
			w.writeAttribute("count", String.valueOf(transcripts.size()));
			for(GTFLine transcript: transcripts) {
				final String transcript_id = transcript.getTranscriptId();
				final CountRows countRow =  transcript2y.get(transcript_id);
				writeStartRecord(w, transcript, countRow, interval);
				
				final List<GTFLine> components = transcript_components.getOrDefault(transcript_id,Collections.emptyList());
				for(GTFLine c: components) {
					write(w,c, countRow, interval);
					}
				
				writeEndRecord(w, transcript);
				}
			w.writeEndElement();//transcripts
			writeEndRecord(w, gene);
			}
		
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		this.skip_attributes_set.addAll(
				Arrays.stream(this.skipAttributeStr.split("[ ;,]+")).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toSet())
				);
		
		XMLStreamWriter w=null;
		PrintWriter fw=null;
		try {
			final String inputName=oneFileOrNull(args);
			final GTFCodec codec = new GTFCodec();
			final SAMSequenceDictionary dict;
			
			if(!StringUtil.isBlank(this.intervalStr) && inputName==null) {
				LOG.info("--interval require tabix indexed file.");
				return -1;
				}
			
			if(this.faidxPath!=null) {
				dict = SequenceDictionaryUtils.extractRequired(this.faidxPath);
				}
			else
				{
				dict=null;
				}
			
			final XMLOutputFactory xof=XMLOutputFactory.newFactory();
			if(this.outputFile==null)
				{
				w = xof.createXMLStreamWriter(stdout(), "UTF-8");
				}
			else
				{
				w = xof.createXMLStreamWriter((fw=super.openPathOrStdoutAsPrintWriter(outputFile)));
				}
			w.writeStartDocument("UTF-8","1.0");
			w.writeStartElement("gtf");
			if(!StringUtil.isBlank(inputName)) {
				w.writeAttribute("filename", inputName);
				}
			final UnaryOperator<String> contigConverter;
			if(dict!=null) {
				new DictionaryXmlSerializer().writeDictionary(w, dict);
				contigConverter = ContigNameConverter.fromOneDictionary(dict);
				codec.setContigNameConverter(contigConverter);
				}
			else {
				contigConverter = ContigNameConverter.getIdentity();
				}
			
			
			
			if(!StringUtil.isBlank(this.intervalStr)) {
				if(inputName==null) {
					LOG.info("--interval require tabix indexed file.");
					return -1;
					}
				
				try(TabixReader tabixR = new TabixReader(inputName)) {
					final SimpleInterval interval = new IntervalParser()
							.contigConverter(contigConverter)
							.apply(this.intervalStr)
							.orElseThrow();
					final List<GTFLine> records = new ArrayList<>();
					final Set<String> gene_ids = new HashSet<>();
					// first pass, look for the min/max of genes
					TabixReader.Iterator iter1 = tabixR.query(interval.getContig(), Math.max(1,interval.getStart()-1), interval.getEnd()+1);
					final MinMaxInteger minmax = new MinMaxInteger(interval.getStart(), interval.getEnd());
					for(;;) {
						final String line = iter1.next();
						if(line==null) break;
						final GTFLine rec = codec.decode(line);
						if(rec==null) continue;
						if(this.disable_group_by_gene) {
							write(w,rec, null,interval);
							}
						else
							{
							if(!rec.isGene() || StringUtils.isBlank(rec.getGeneId())) continue;
							if(!rec.overlaps(interval)) continue;// the interval was extended...
							minmax.accept(rec.getStart());
							minmax.accept(rec.getEnd());
							gene_ids.add(rec.getGeneId());
							}
						}
					if(!minmax.isEmpty()) {
						iter1 = tabixR.query(interval.getContig(), Math.max(1,minmax.getMinAsInt()-1), minmax.getMaxAsInt()+1);
						for(;;) {
							final String line = iter1.next();
							if(line==null) break;
							final GTFLine rec = codec.decode(line);
							final String gene_id = rec.getGeneId();
							if(StringUtils.isBlank(gene_id)) continue;
							if(!gene_ids.contains(gene_id)) continue;
							records.add(rec);
							}
						w.writeStartElement("contig");
						w.writeAttribute("name", interval.getContig());
						w.writeAttribute("start",String.valueOf(interval.getStart()));
						w.writeAttribute("end",String.valueOf(interval.getEnd()));
						compile(w, records, interval);
						w.writeEndElement();//contig
						}
					}
				}
			else if(this.disable_group_by_gene) {
				try(BufferedReader br=StringUtil.isBlank(inputName)?
						IOUtils.openStreamForBufferedReader(stdin()):
						IOUtils.openURIForBufferedReading(inputName)
						) {
					for(;;) {
						final String line = br.readLine();
						if(line==null) break;
						if(line.startsWith("#") || StringUtil.isBlank(line)) continue;
						final GTFLine rec = line==null?null:codec.decode(line);
						if(rec==null) continue;
						write(w,rec, null, null);
						}
					}
				}
			else
				{
				final OrderChecker<GTFLine> orderChecker = new OrderChecker<>();
				try(BufferedReader br=StringUtil.isBlank(inputName)?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(inputName)
					) {
					final List<GTFLine> records = new ArrayList<>();
					String prevContig=null;
					for(;;) {
						final String line = br.readLine();
						if(line==null) break;
						if(line.startsWith("#") || StringUtil.isBlank(line)) continue;
						final GTFLine rec = codec.decode(line);
						if(rec==null) continue;
						orderChecker.apply(rec);
						if(prevContig!=null && !prevContig.equals(rec.getContig())) {
							if(!records.isEmpty()) {
								w.writeStartElement("contig");
								w.writeAttribute("name", prevContig);
								compile(w, records, null);
								records.clear();								
								w.writeEndElement();//contig
								}
							}
						records.add(rec);
						prevContig = rec.getContig();
						}
					
					if(!records.isEmpty()) {
						w.writeStartElement("contig");
						w.writeAttribute("name", prevContig);
						compile(w, records, null);
						w.writeEndElement();//contig
						}
					}
				}
			
			w.writeEndElement();
			w.writeEndDocument();
			w.flush();
			return 0;
			}
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			if(fw!=null) fw.close();
			}
		}

	
	public static void main(final String[] args) throws IOException
		{
		new Gtf2Xml().instanceMainWithExit(args);
		}
}
