/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import gov.nih.nlm.ncbi.gb.GBFeature;
import gov.nih.nlm.ncbi.gb.GBFeatureQuals;
import gov.nih.nlm.ncbi.gb.GBQualifier;
import gov.nih.nlm.ncbi.gb.GBSet;
import gov.nih.nlm.ncbi.gb.ObjectFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
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
description="Experimental genbank to GFF",
keywords={"xml","ncbi","genbank","convert","gff","gb"}
)
public class GenbankToGff3 extends Launcher {
	private static final Logger LOG = Logger.build(GenbankToGff3.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@SuppressWarnings("unused")
	private static final ObjectFactory _fool_javac = null;
	private final Set<String> headers = new HashSet<>();
	private PrintWriter tmpWriter = null;
	private final Unmarshaller unmarshaller;
	private long ID_GENERATOR = 0L;
	private final List<MyFeature> buffer= new ArrayList<>();
	
	/** wrapper around GBFeature */
	private class MyFeature
		{
		private final long id;
		private final GBFeature feat;
		MyFeature(final GBFeature feat)
			{
			this.feat = feat;
			this.id = (++ID_GENERATOR);
			//prevent NPE
			if(this.feat.getGBFeatureQuals()==null)
				{	
				this.feat.setGBFeatureQuals(new GBFeatureQuals());
				}
			
			this.feat.getGBFeatureQuals().
				getGBQualifier().
				removeIf(Q->(Q.getGBQualifierName().equals("translation") || 
						     Q.getGBQualifierName().equals("transcription") ));
		
			}
		public String getType() 
			{
			return this.feat.getGBFeatureKey().
					replace("5'", "five_prime_").
					replace("3'", "three_prime_")
					;
			}
		public Map<String,String> getQualMap() {
			
			final Map<String,String> hash = new HashMap<String,String>();
			// do not use java stream because there are some dups'
			for(final GBQualifier G:feat. getGBFeatureQuals().getGBQualifier())
				{
				hash.put(G.getGBQualifierName(),G.getGBQualifierValue());
				}
			return hash;
			}
		
		private String findQualifier(final String key)
			{
			return this.feat.getGBFeatureQuals().getGBQualifier().stream().
					filter(G->G.getGBQualifierName().equals(key)).
					map(G->G.getGBQualifierValue()).
					findFirst().
					orElse(null);
			}
		
		public void print() {
			final Map<String,String> qualifiers= getQualMap();
			
			if(getType().equals("source")) {
				final String taxon = findQualifier("db_xref");
				if(taxon!=null && taxon.startsWith("taxon:"))
					{
					GenbankToGff3.this.headers.add("species https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id="+taxon.substring(6).trim());
					}
				getIntervals().stream().forEach(I-> {
					GenbankToGff3.this.headers.add("sequence-region "+I.getContig()+" "+I.getStart()+" "+I.getEnd());	
					
					});
				}
			
				getIntervals().stream().
					forEach(I->print(I,qualifiers));
				}
			
			public String getId() {
				return getType()+this.id;
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

			
			void print(
					final Interval interval,
					Map<String,String> quals
					)
				{
				//reset ID and parentID
				quals.remove("ID");
				quals.remove("ParentID");
				
				Integer phase = null;
				if( getType().equals("CDS")) {
					final String q= this.findQualifier("codon_start");
					if(q!=null && q.matches("[0-9]+"))
						{
						phase= Integer.parseInt(q);
						phase=phase-1;
						}
					}
				final String locus_tag = this.findQualifier("locus_tag");
				if(!StringUtil.isBlank(locus_tag))
					{
					quals.put("ID", this.getId());
					
					if(!getType().equals("gene")) {
						final String parentId=GenbankToGff3.this.buffer.
							stream().
							filter(G->G.getType().equals("gene") && locus_tag.equals(G.findQualifier("locus_tag"))).
							map(G->G.getId()).
							findFirst().
							orElse(null);
						if(!StringUtil.isBlank(parentId))
							{
							quals.put("ParentID",parentId);
							}
						}
					}
				
				final PrintWriter w  = GenbankToGff3.this.tmpWriter;
				w.print(interval.getContig());
				w.print('\t');
				w.print("genbank");//source
				w.print('\t');
				w.print(getType());
				w.print('\t');
				w.print(interval.getStart());
				w.print('\t');
				w.print(interval.getEnd());
				w.print('\t');
				w.print('.');//score
				w.print('\t');
				w.print(interval.isNegativeStrand()?'-':'+');
				w.print('\t');
				w.print(phase==null?".":String.valueOf(phase));
				w.print('\t');
				final String qualifiers = quals.keySet().stream().
						map(K->K+"="+escape(quals.get(K))).
						collect(Collectors.joining(";"))
						;
				
				
				if(StringUtil.isBlank(qualifiers))
					{
					w.print('.');
					}
				else
					{
					w.print(qualifiers);
					}
				
				w.println();
				}
		
		public List<Interval> getIntervals()
			{
			final List<Interval> intervals = new ArrayList<>();
			this.feat.getGBFeatureIntervals().getGBInterval().forEach(T->
				{
				final String acn = T.getGBIntervalAccession();
				if(StringUtil.isBlank(acn))
					{
					LOG.warn("interval acn missing");
					return;
					}
				
				final int intervalFrom ;
				final int intervalTo ;
				
				if(!StringUtil.isBlank(T.getGBIntervalPoint()))
					{
					intervalFrom = Integer.parseInt(T.getGBIntervalPoint());
					intervalTo=intervalFrom;
					}
					
				else
					{
					intervalFrom = Integer.parseInt(T.getGBIntervalFrom());
					intervalTo = Integer.parseInt(T.getGBIntervalTo());
					}
				final Interval interval =  new Interval(
						acn,
						Math.min(intervalFrom,intervalTo),
						Math.max(intervalFrom,intervalTo),
						intervalFrom>intervalTo,
						"gbfeature"
						);
				intervals.add(interval);
				});
			
			return intervals;
			}
		}
	
	public GenbankToGff3() {
		
		
		JAXBContext jc;
		try {
			jc = JAXBContext.newInstance(GBSet.class.getPackage().getName());
			this.unmarshaller = jc.createUnmarshaller();
		} catch (final JAXBException e) {
			LOG.error(e);
			throw new RuntimeException(e);
			}
		}
	
	
	
	

	
	private void dump() {
		this.buffer.forEach(B->B.print());
		this.buffer.clear();
		}	
	private void parseGenBank(final XMLEventReader r) throws JAXBException,XMLStreamException,IOException{
		
		while(r.hasNext())
			{
			final XMLEvent evt= r.peek();
			if(evt.isStartElement())
				{
				final String name = evt.asStartElement().getName().getLocalPart();
				 if(name.equals("GBFeature"))
					{
					final GBFeature feature = this.unmarshaller.unmarshal(r, GBFeature.class).getValue(); 
					final MyFeature my = new MyFeature(feature);
					if(!my.getIntervals().isEmpty()) {
						// we can print things like SNP right now...
						if(StringUtil.isBlank( my.findQualifier("locus_tag")))
							{
							my.print();
							}
						else
							{
							this.buffer.add(my);			
							}
						}
					continue;
					}
				 else if(name.equals("GBSeq"))
					{
					r.nextEvent();
					dump();
					continue;
					}
				}
			//consumme event
			r.nextEvent();
			}
		dump();
		}
	
	@Override
	public int doWork(final List<String> args) {
		final List<Path> inputs = IOUtil.unrollFiles(
				args.stream().map(S->new File(S)).collect(Collectors.toSet()),
				".xml",".gb",".genbank"
				).stream().map(F->F.toPath()).collect(Collectors.toList());
		Path tmpPathBody= null;
		BufferedWriter bw = null;
		PrintWriter w=null;
		final XMLInputFactory xif = XMLInputFactory.newInstance();
		try {
			
			tmpPathBody = Files.createTempFile("tmp.", ".gff3");
			this.tmpWriter = new PrintWriter( Files.newBufferedWriter(tmpPathBody));
			
			if(inputs.isEmpty()) {
				LOG.info("reading stdin");
				final Reader r = new InputStreamReader(stdin());
				final XMLEventReader xr = xif.createXMLEventReader(r);
				parseGenBank(xr);
				xr.close();
				r.close();
				} 
			else for(final Path path:inputs)
				{
				LOG.info("reading "+path);
				final Reader r = Files.newBufferedReader(path);
				final XMLEventReader xr = xif.createXMLEventReader(r);
				parseGenBank(xr);
				xr.close();
				r.close();
				}
			
			
			this.tmpWriter.flush();
			this.tmpWriter.close();
			this.tmpWriter=null;
			
			w = super.openFileOrStdoutAsPrintWriter(outputFile);
			for(final String h:headers) w.println("##"+h);
			
			IOUtils.copyTo(tmpPathBody, w);
			w.flush();
			w.close();
			w=null;
			return 0;
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(bw);
			CloserUtil.close(this.tmpWriter);
			
			if(tmpPathBody!=null) try { Files.deleteIfExists(tmpPathBody);}
				catch(Exception err) {}
			headers.clear();
			}
		}
	
public static void main(final String[] args) {
	new GenbankToGff3().instanceMainWithExit(args);
	}
}
