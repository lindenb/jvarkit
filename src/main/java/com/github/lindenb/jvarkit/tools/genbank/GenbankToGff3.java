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
package com.github.lindenb.jvarkit.tools.genbank;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.annotation.Strand;
import htsjdk.tribble.gff.Gff3BaseData;
import htsjdk.tribble.gff.Gff3Constants;
/**
BEGIN_DOC

## Warnings

 - still experimental
 - efetch doesn't always work: https://gist.github.com/lindenb/6c4f36bdf29a3108e103e3a5a0b1aff7

## Example

```
$ wget -O - -q  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=AF338247.1,D38149.1&rettype=gbwithparts&retmode=xml&api_key=62b713e0cd85e6ac79699ecdfa72e85af009"  | java -jar dist/gb2gff.jar 

##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10970
##species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=10941
##sequence-region AF338247.1 1 2032
##sequence-region D38149.1 1 1087
AF338247.1	genbank	source	1	2032	.	+	.	strain=M;db_xref=taxon:10941;note="rearranged segment 7";organism="Human rotavirus A";segment=7R;clone=M1;mol_type="genomic RNA"
AF338247.1	genbank	five_prime_UTR	1	34	.	+	.	.
AF338247.1	genbank	CDS	35	967	.	+	1	codon_start=1;product=NSP3;protein_id=AAK74117.1;transl_table=1
AF338247.1	genbank	CDS	993	1925	.	+	1	note="duplicated ORF";codon_start=1;product=NSP3;protein_id=AAK74118.1;transl_table=1
AF338247.1	genbank	three_prime_UTR	1926	2032	.	+	.	.
D38149.1	genbank	source	1	1087	.	+	.	strain=A5-16;db_xref=taxon:10970;organism="Rotavirus sp.";segment=5;mol_type="genomic RNA"
D38149.1	genbank	gene	33	185	.	+	.	gene=NSP1
D38149.1	genbank	CDS	33	185	.	+	1	codon_start=1;protein_id=BAA07347.1;gene=NSP1;transl_table=1
D38149.1	genbank	variation	141	141	.	+	.	note="bases 142-641 in D38148 deleted";gene=NSP1
```

## input 

Input is a list of *.xml genbank files (DTD/Schema: https://www.ncbi.nlm.nih.gov/dtd/NCBI_GBSeq.dtd ). Otherwise assume that file is a list of filenames and unfold it.

END_DOC
 */
@Program(name="gb2gff",
description="Experimental Genbank XML to GFF",
keywords={"xml","ncbi","genbank","convert","gff","gb"},
modificationDate="20210429",
creationDate="20180215",
generate_doc = false
)
public class GenbankToGff3 extends Launcher {
	private static final Logger LOG = Logger.build(GenbankToGff3.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	private static long ID_GENERATOR = 0L;
	private final static double NO_SCORE=-1;
	private final static int NO_PHASE=-1;
	private final static String DUMMY_CTG="chr_";

	
	private static class GBFeature implements Locatable {
		String key;
		final List<Interval> gbIntervals = new ArrayList<>();
		final Map<String,String> qualifiers= new HashMap<>();
		public boolean isGene() {
			return key.equals("gene");
			}
		public boolean isTranscriptOf(final String g)  {
				return isTranscript() && qualifiers.getOrDefault("gene","").equals(g);
				}
		public boolean isTranscript() { 
			return key.equals("misc_RNA") || key.equals("mRNA")|| key.equals("ncRNA");
		}
		public boolean isExon()  { return key.equals("exon");}
		public boolean isCDS()  { return key.equals("CDS");}

		
		public String getGeneId() {
			if(qualifiers.containsKey("gene")) {
				return qualifiers.get("gene");
				}
			if(qualifiers.containsKey("locus_tag")) {
				return qualifiers.get("locus_tag");
				}
			return null;
			}
		@Override
		public String getContig() {
			return DUMMY_CTG;
			}
		@Override
		public int getStart() {
			return gbIntervals.stream().mapToInt(R->R.getStart()).min().getAsInt();
			}
		@Override
		public int getEnd() {
			return gbIntervals.stream().mapToInt(R->R.getEnd()).max().getAsInt();
			}
		public Strand getStrand() {
			return gbIntervals.stream().anyMatch(R->R.isNegativeStrand())?Strand.NEGATIVE:Strand.FORWARD;
			}
		}
	
	
	private class Entry {
		String name = "undefined";
		int length=-1;
		final List<GBFeature> features = new ArrayList<>();
		}
	

	
	
	private String escape(final String s)
		{
		if(s.contains(" ") || s.contains("=") || s.contains("\\")|| s.contains(";") || s.contains("\""))
			{
			final StringBuilder sb=new StringBuilder("\"");
			for(int i=0;i< s.length();++i)
				{
				switch(s.charAt(i))
					{
					case '\\': sb.append("\\\\"); break;
					case '\'': sb.append("\\\'"); break;
					case '\"': sb.append("\\\""); break;
					default: sb.append(s.charAt(i));break;
					}
				}
			sb.append("\"");
			return sb.toString();
			}
		return s;
		}
	
	
	

	
	private Interval parseGBInterval(final XMLEventReader r) throws XMLStreamException,IOException{
		int start=-1;
		int end=1;
		String accession="";
		while(r.hasNext())
			{
			final XMLEvent evt= r.nextEvent();
			if(evt.isStartElement()) {
				final String name = evt.asStartElement().getName().getLocalPart();
				if(name.equals("GBInterval_from")) {
					start = Integer.parseInt(r.getElementText());
					}
				else if(name.equals("GBInterval_to")) {
					end = Integer.parseInt(r.getElementText());
					}
				else if(name.equals("GBInterval_point")) {
					start = Integer.parseInt(r.getElementText());
					end = start;
					}
				else if(name.equals("GBInterval_accession")) {
					accession = r.getElementText();
					}
				}
			else if(evt.isEndElement()) {
				final String name = evt.asEndElement().getName().getLocalPart();
				if(name.equals("GBInterval")) break;
				}
			}
		
		return new Interval(DUMMY_CTG,Math.min(start, end),Math.max(start, end),start>end,accession);
		}
	private Map.Entry<String, String> parseGBQualifier(final XMLEventReader r) throws XMLStreamException,IOException{
		String key="";
		String value="true";
		while(r.hasNext())
			{
			final XMLEvent evt= r.nextEvent();
			if(evt.isStartElement()) {
				final String name = evt.asStartElement().getName().getLocalPart();
				if(name.equals("GBQualifier_name")) {
					key = r.getElementText();
					}
				else if(name.equals("GBQualifier_value")) {
					value = r.getElementText();
					}
				}
			else if(evt.isEndElement()) {
				final String name = evt.asEndElement().getName().getLocalPart();
				if(name.equals("GBQualifier")) break;
				}
			}
		return new AbstractMap.SimpleEntry<>(key,value);
		}

	
	private GBFeature parseFeature(final XMLEventReader r) throws XMLStreamException,IOException{
		final GBFeature feature = new GBFeature();
		while(r.hasNext())
			{
			final XMLEvent evt= r.nextEvent();
			if(evt.isStartElement()) {
				final String name = evt.asStartElement().getName().getLocalPart();
				if(name.equals("GBFeature_key")) {
					feature.key = r.getElementText().
							replace("5'", "five_prime_").
							replace("3'", "three_prime_");
					}
				else if(name.equals("GBInterval")) {
					feature.gbIntervals.add(parseGBInterval(r));
					}
				else if(name.equals("GBQualifier")) {
					final Map.Entry<String, String> kv = parseGBQualifier(r);
					if(!(kv.getKey().equals("transcription") || kv.getKey().equals("translation") || kv.getKey().equals("note") || kv.getKey().equals("peptide"))) {
						feature.qualifiers.put(kv.getKey(),kv.getValue());
						}
					}
				}
			else if(evt.isEndElement()) {
				final String name = evt.asEndElement().getName().getLocalPart();
				if(name.equals("GBFeature")) {
					break;
					}
				}
			}
		
		return feature;
		}
	
	private List<GBFeature> parseFeatureTable(final XMLEventReader r) throws XMLStreamException,IOException {
		final List<GBFeature> L = new ArrayList<>();
		while(r.hasNext())
			{
			final XMLEvent evt= r.nextEvent();
			if(evt.isStartElement()) {
				final String name = evt.asStartElement().getName().getLocalPart();
				if(name.equals("GBFeature")) {
					final GBFeature feature =  parseFeature(r);
					L.add(feature);
					}
				}
			else if(evt.isEndElement()) {
				final String name = evt.asEndElement().getName().getLocalPart();
				if(name.equals("GBSeq_feature-table")) break;
				}
			}
		Collections.sort(L,(A,B)->Integer.compare(A.getStart(), B.getStart()));
		return L;
		}

	
	private Entry parseGBSeq(final XMLEventReader r) throws XMLStreamException,IOException{
		final Entry entry = new Entry();
		while(r.hasNext())
			{
			final XMLEvent evt= r.nextEvent();
			if(evt.isStartElement()) {
				final String name = evt.asStartElement().getName().getLocalPart();
				if(name.equals("GBSeq_locus")) {
					entry.name  = r.getElementText();
					}
				else if(name.equals("GBSeq_length")) {
					entry.length  = Integer.parseInt(r.getElementText());
					}
				else if(name.equals("GBSeq_feature-table")) {
					entry.features.addAll(parseFeatureTable(r));
					}
				
			}
			else if(evt.isEndElement()) {
				final String name = evt.asEndElement().getName().getLocalPart();
				if(name.equals("GBFeature")) break;
				}
			}
		return entry;
		}
	
	private final List<Entry> parseGenBank(final XMLEventReader r) throws XMLStreamException,IOException{
		final List<Entry> L = new ArrayList<>();
		while(r.hasNext())
			{
			final XMLEvent evt= r.nextEvent();
			if(evt.isStartElement()) {
				final String name = evt.asStartElement().getName().getLocalPart();
				if(name.equals("GBSeq")) {
					L.add(this.parseGBSeq(r));
					}
				}
			}
		return L;
		}
	private void print(PrintWriter pw,final Gff3BaseData record) {
		pw.print(record.getContig());
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(record.getSource());
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(record.getType());
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(record.getStart());
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(record.getEnd());
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(Gff3Constants.UNDEFINED_FIELD_VALUE);//score
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(record.getStrand().encodeAsChar());
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(record.getPhase()==NO_PHASE?Gff3Constants.UNDEFINED_FIELD_VALUE:String.valueOf(record.getPhase()));
		pw.print(Gff3Constants.FIELD_DELIMITER);
		pw.print(
				record.getAttributes().entrySet().stream().
				map(KV->KV.getKey()+Gff3Constants.KEY_VALUE_SEPARATOR+escape(KV.getValue().get(0))).collect(Collectors.joining(""+Gff3Constants.ATTRIBUTE_DELIMITER)));
		pw.println();
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<Path> inputs = IOUtil.unrollFiles(
				args.stream().map(S->new File(S)).collect(Collectors.toSet()),
				".xml",".gb",".genbank"
				).stream().map(F->F.toPath()).collect(Collectors.toList());
		
		final XMLInputFactory xif = XMLInputFactory.newInstance();
		try {
			
			final List<Entry> L = new ArrayList<>();
			if(inputs.isEmpty()) {
				LOG.info("reading stdin");
				try(final Reader r = new InputStreamReader(stdin())) {
					final XMLEventReader xr = xif.createXMLEventReader(r);
					L.addAll(parseGenBank(xr));
					xr.close();
					}
				} 
			else for(final Path path:inputs)
				{
				LOG.info("reading "+path);
				try(final Reader r = Files.newBufferedReader(path)) {
					final XMLEventReader xr = xif.createXMLEventReader(r);
					L.addAll(parseGenBank(xr));
					xr.close();
					}
				}
			
			

			try(PrintWriter w = super.openFileOrStdoutAsPrintWriter(outputFile)) {
				
				w.println("##gff-version 3");
				w.println("##format: gff3");
				for(Entry entry: L) {
					w.println("##sequence-region "+entry.name+" "+entry.length);
					}
				for(Entry entry: L) {
					Map<String,List<String>> attCtg = Collections.singletonMap("ID", Collections.singletonList(entry.name));
					Gff3BaseData gffCtg = new Gff3BaseData(
							entry.name,
							"GB",
							"chromosome",
							1,
							entry.length,
							NO_SCORE,
							Strand.NONE,
							NO_PHASE,
							attCtg
							);
					print(w,gffCtg);
										
							
					for(GBFeature gene: entry.features.stream().filter(G->G.isGene()).collect(Collectors.toList())) {
						final Map<String,List<String>> gattributes = new HashMap<>();
						String gene_id = "G"+(++ID_GENERATOR);
						final String geneName = gene.qualifiers.getOrDefault("gene", "");
						if(StringUtils.isBlank(geneName)) {
							LOG.warn("no gene name for "+gene);
							continue;
							}
						gene.qualifiers.entrySet().forEach(KV->gattributes.put(KV.getKey(),Collections.singletonList(KV.getValue())));
						
						gattributes.put("Parent",Collections.singletonList(entry.name));
						gattributes.put("ID",Collections.singletonList(gene_id));
						gattributes.put("Name",Collections.singletonList(geneName));
						
						
						final Gff3BaseData gffGene = new Gff3BaseData(
								entry.name,
								"GB",
								gene.key,
								gene.getStart(),
								gene.getEnd(),
								NO_SCORE,
								gene.getStrand(),
								NO_PHASE,
								gattributes
								);
						print(w,gffGene);
						
						for(GBFeature transcript: entry.features.stream().filter(G->G.isTranscriptOf(geneName)).collect(Collectors.toList())) {
							final Map<String,List<String>> tattributes = new HashMap<>();
							String transcript_id = "T"+(++ID_GENERATOR);
							transcript.qualifiers.entrySet().forEach(KV->tattributes.put(KV.getKey(),Collections.singletonList(KV.getValue())));
							
							tattributes.put("Parent",Collections.singletonList(gene_id));
							tattributes.put("ID",Collections.singletonList(transcript_id));
							//tattributes.put("Name",Collections.singletonList(geneName));
							
							
							final Gff3BaseData gffRNA = new Gff3BaseData(
									entry.name,
									"GB",
									transcript.key,
									transcript.getStart(),
									transcript.getEnd(),
									NO_SCORE,
									transcript.getStrand(),
									NO_PHASE,
									tattributes
									);
							print(w,gffRNA);
							}
						
						}
					}
				w.flush();
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
public static void main(final String[] args) {
	new GenbankToGff3().instanceMainWithExit(args);
	}
}
