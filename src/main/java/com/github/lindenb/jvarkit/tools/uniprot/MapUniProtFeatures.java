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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.uniprot;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.StringWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptReader;
import com.github.lindenb.jvarkit.uniprot.Uniprot;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;

/**

BEGIN_DOC

## Warning

this program is not fully tested. Please check the results

##Example

```bash

wget -O knownGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV47.txt.gz
wget -O knownGene.sql http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/wgEncodeGencodeBasicV47.sql


$ java  -jar dist/jvarkit mapuniprot \
	-R /path/to/hg38.fasta \
	-u /path/uri/uniprot.org/uniprot_sprot.xml.gz  \
	-k knownGene.txt.gz | gunzip -c | awk -F '        ' '{if($2 ~ ".*_.*") next; OFS="       "; gsub(/chr/,"",$2);print;}'   ) |\
	LC_ALL=C sort -t '	' -k1,1 -k2,2n -k3,3n  | uniq | head


1	69090	69144	topological_domain	1000	+	69090	69144	255,0,0	1	54	0
1	69144	69216	transmembrane_region	1000	+	69144	69216	255,0,0	1	72	0
1	69216	69240	topological_domain	1000	+	69216	69240	255,0,0	1	24	0
1	69240	69306	transmembrane_region	1000	+	69240	69306	255,0,0	1	66	0
1	69306	69369	topological_domain	1000	+	69306	69369	255,0,0	1	63	0
1	69357	69636	disulfide_bond	1000	+	69357	69636	255,0,0	1	279	0
1	69369	69429	transmembrane_region	1000	+	69369	69429	255,0,0	1	60	0
1	69429	69486	topological_domain	1000	+	69429	69486	255,0,0	1	57	0
1	69486	69543	transmembrane_region	1000	+	69486	69543	255,0,0	1	57	0
1	69543	69654	topological_domain	1000	+	69543	69654	255,0,0	1	111	0
```
END_DOC
 
 */
@Program(
	name="mapuniprot",
	description="map uniprot features on reference genome.",
	keywords={"uniprot","bed","fasta","reference","xml"},
	modificationDate = "20250331",
	jvarkit_amalgamion = true
	)
public class MapUniProtFeatures extends Launcher
	{
	private static Logger LOG=Logger.of(MapUniProtFeatures.class);
    private static class Range
    	{
    	Range(final int start,final int end)
    		{
    		this.start=start;
    		this.end=end;
    		}
    	int start;
    	int end;
    	}
	private  class UBed
		{
		String contig;
		String featureType;
		String featureDesc;
		List<Range> positions=new ArrayList<Range>();
		byte strand;
		String name;
		
		String getContig()
			{
			return contig;
			}
		
		int start()
			{
			return this.positions.get(0).start;
			}
		int end()
			{
			return this.positions.get(positions.size()-1).end;
			}
		@Override
		public String toString() {
			return toBed12();
			}

		
		private String toBed12()
			{
			final StringBuilder b=new StringBuilder();
			b.append(getContig());
			b.append('\t').append(start());
			b.append('\t').append(end());
			b.append('\t').append(name.replaceAll("[ ]+", "_"));
			b.append('\t').append(1000);
			b.append('\t').append((char)strand);
			b.append('\t').append(start());
			b.append('\t').append(end());
			b.append('\t').append("255,0,0");
			b.append('\t').append(positions.size());
			b.append('\t');
			
			for(int i=0;i< positions.size();++i)
				{
				if(i>0) b.append(',');
				b.append(positions.get(i).end-positions.get(i).start);
				}
			b.append('\t');
			for(int i=0;i< positions.size();++i)
				{
				if(i>0) b.append(',');
				b.append(positions.get(i).start-start());
				}
			b.append('\t');
			b.append(this.featureType);
			return b.toString();
			}
		
		private void toAnnotate(PrintWriter pw) {
			for(Range r: this.positions) {
				pw.print(this.getContig());
				pw.print('\t');
				pw.print(r.start);
				pw.print('\t');
				pw.print(r.end);
				pw.print('\t');
				pw.print(name.replaceAll("[ ]+", "_"));
				pw.print('\t');
				String s= featureType.replaceAll("[ ;=]+", "_");
				if(StringUtils.isBlank(s)) s=".";
				pw.print(s);
				pw.print('\t');
				s= featureDesc.replaceAll("[ ;=]+", "_");
				if(StringUtils.isBlank(s)) s=".";
				pw.print(s);
				pw.println();
				}
			}
		
		void print(FormatOut fmt,PrintWriter pw)
			{
			switch(fmt) {
				case annotate:
					this.toAnnotate(pw);
					break;
				case bed12: //
				default:
					pw.println(this.toBed12());
					break;
				}
			}
		}
	
	private enum FormatOut { annotate, bed12};
	
	@Parameter(names={"-k","--transcripts","--genpred"},description=UcscTranscriptReader.OPT_DESC,required=true)
	private String kgUri=null;
	
	@Parameter(names="-R",description=Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION, required = true)
	private Path referencePath = null;
	
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile =null;

	@Parameter(names="--format",description="output format 'annotate' is suitable for bcftools annotate")
	private FormatOut formatOut =FormatOut.bed12;

	@Parameter(names="--debug",description="debug",hidden = true)
	private boolean debug=false;

	private UcscTranscript.Peptide getPeptide(final UcscTranscript kg,final GenomicSequence genomic) {
		UcscTranscript.MessengerRNA mRNA= kg.getMessengerRNA(genomic);
		UcscTranscript.CodingRNA cDNA= mRNA.getCodingRNA();
		return cDNA.getPeptide();
		}
	
	private void recursiveBuildXml(final XMLEventReader xr,final Document dom,final Node root, boolean append,int depth) throws XMLStreamException{
		while(xr.hasNext()) {
			XMLEvent evt = xr.nextEvent();
			switch(evt.getEventType()) {
				case XMLEvent.START_ELEMENT:
					final StartElement E = evt.asStartElement();
					final String lclName = E.getName().getLocalPart();
					
					final Element n;
					boolean appen_child = append;
					if(append) {
						if(lclName.equals("reference")) appen_child = false;
						else if(lclName.equals("dbReference")) appen_child = false;
						else if(lclName.equals("proteinExistence")) appen_child = false;
						else if(lclName.equals("keyword")) appen_child = false;
						else if(lclName.equals("evidence")) appen_child = false;
						}
					
					if(appen_child) {
						n = dom.createElement(lclName);
						root.appendChild(n);
						
						for(Iterator<?> iter =E.getAttributes();iter.hasNext();) {
							Attribute att =(Attribute)iter.next();
							n.setAttribute(att.getName().getLocalPart(), att.getValue());
							}
						}
					else
						{
						n = null;
						}
					recursiveBuildXml(xr,dom,n, appen_child,depth+1);
					break;
				case XMLEvent.END_ELEMENT:
					return;
				case XMLEvent.END_DOCUMENT:
					throw new IllegalStateException();
				case XMLEvent.CDATA:
				case XMLEvent.CHARACTERS:
					if(append) root.appendChild(dom.createTextNode(evt.asCharacters().getData()));
					break;
				default:
					break;
				}
		}
	}
	
	
	private String toString(Document dom) {
		try {
			StringWriter  sw = new StringWriter();
			Transformer tr= TransformerFactory.newInstance().newTransformer();
			tr.transform(new DOMSource(dom), new StreamResult(sw));
			return sw.toString();
		} catch(Throwable err) {
			return "#ERROR";
		}
	}
	
	@Override
	public int doWork(List<String> args)
		{
		final Map<String,List<UcscTranscript>> protid2transcripts = new HashMap<>(100_000);

		try(ReferenceSequenceFile indexedFastaSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referencePath) )
			{
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(indexedFastaSequenceFile);
			String line;
			final UcscTranscriptCodec codec = new UcscTranscriptCodec(kgUri);
			try(BufferedReader r=IOUtils.openURIForBufferedReading(kgUri)) {
				while((line=r.readLine())!=null)
					{
					final UcscTranscript kg = codec.decode(line);
					if(kg==null || !kg.isProteinCoding()) continue;
					if(dict.getSequence(kg.getContig())==null || StringUtils.isBlank(kg.getName2())) continue;
					
					List<UcscTranscript> L=protid2transcripts.get(kg.getName2());
					if(L==null)
						{
						L=new ArrayList<>();
						protid2transcripts.put(kg.getName2(), L);
						}
					
					L.add(kg);
					}
				}
			
			
			
			final XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
			xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
			xmlInputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.FALSE);
			
			xmlInputFactory.setXMLResolver(new XMLResolver() {
					@Override
					public Object resolveEntity(String arg0, String arg1, String arg2,
							String arg3) throws XMLStreamException
						{
						LOG.info("resolveEntity:" +arg0+"/"+arg1+"/"+arg2);
						return null;
					}
				});
			
			//SortingCollection<UBed> mappedFeatures=SortingCollection.newInstance(UBed.class, new UBedCodec(),new UBedCmp(),super.MAX_RECORDS_IN_RAM);
			//mappedFeatures.setDestructiveIteration(true);
			
			
			try(Reader fr = super.openBufferedReader(this.oneFileOrNull(args))) {
				final XMLEventReader rx=xmlInputFactory.createXMLEventReader(fr);
				final DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
				dbf.setNamespaceAware(false);
				final DocumentBuilder db = dbf.newDocumentBuilder();
				
				final XPath xpath = XPathFactory.newInstance().newXPath();
				
				try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
					GenomicSequence genomic=null;
					while(rx.hasNext())
						{
						final XMLEvent evt=rx.nextEvent();
						if(!(evt.isStartElement() && Uniprot.isNamespace(evt.asStartElement().getName().getNamespaceURI()) && evt.asStartElement().getName().getLocalPart().equals("entry")))
							{
							if(debug && evt.isStartElement()) {
								LOG.debug("skipping "+evt.asStartElement().getName());
								}
							continue;
							}
						
						if(debug) {
							LOG.debug("got "+evt.asStartElement().getName());
							}
						
						final Document dom = db.newDocument();
						final Element uniprotNode = dom.createElement("uniprot");
						dom.appendChild(uniprotNode);
						final Element entryNode = dom.createElement("entry");
						uniprotNode.appendChild(entryNode);
						recursiveBuildXml(rx, dom, entryNode,true,0);
						
						
						if(debug) {
							LOG.info("DONE scanning");
							}
						
						if(!((Boolean)xpath.evaluate("/uniprot/entry/feature",dom,XPathConstants.BOOLEAN))) {
							if(debug) {
								LOG.info("no /uniprot/entry/feature found for "+ xpath.evaluate("/uniprot/entry/name",dom,XPathConstants.STRING)+ " "+toString(dom));
								}
							continue;
							}
						
						final Set<String> transcriptIdSet = new HashSet<>();
						NodeList nodeList= (NodeList)xpath.evaluate(
								"/uniprot/entry/accession|"
								+ "/uniprot/entry/name|"
								+ "/uniprot/entry/gene/name|"
								+ "/uniprot/entry/dbReference[type=\"Ensembl\"]/@Id",
								dom,XPathConstants.NODESET);
						for(int i=0;i< nodeList.getLength();i++)
							{
							final String acn = nodeList.item(i).getTextContent();
							if(StringUtils.isBlank(acn) || !protid2transcripts.containsKey(acn)) continue;
							transcriptIdSet.add(acn);
							}
						
						final List<UcscTranscript> genes=transcriptIdSet.stream().
								flatMap(K->protid2transcripts.get(K).stream()).
								collect(Collectors.toList());
						
						if(genes.isEmpty()) {
							if(debug && !transcriptIdSet.isEmpty()) {
								LOG.info("no gene found for "+ String.join(",", transcriptIdSet));
								}
							continue;
							}
						LOG.info("got "+genes+" for "+ xpath.evaluate("/uniprot/entry/name",dom,XPathConstants.STRING));
						final String entrySequence = (String)xpath.evaluate("/uniprot/entry/sequence",dom,XPathConstants.STRING);
						if(StringUtils.isBlank(entrySequence)) {
							if(debug) {
								LOG.info("no sequence found found for "+ xpath.evaluate("/uniprot/entry/name",dom,XPathConstants.STRING));
								}
							continue;
							}
						
						for(UcscTranscript kg:genes)
							{
							if(genomic==null ||  !genomic.hasName(kg.getContig()))
								{
								genomic=new GenomicSequence(indexedFastaSequenceFile, kg.getContig());
								}
							
							if(debug) {
								LOG.debug("testing "+kg.getTranscriptId());
								}
							final UcscTranscript.Peptide pep=getPeptide(kg,genomic);
							if(debug) {
								System.err.println(kg.getTranscriptId()+" peptide: "+pep);
							}
							
							/* scan both sequence while they are the same */
							int sameSequenceLength=0;
							while(  sameSequenceLength < entrySequence.length() &&
									sameSequenceLength < pep.length() 
									)
								{
								if(Character.toUpperCase(entrySequence.charAt(sameSequenceLength))!=Character.toUpperCase(pep.charAt(sameSequenceLength)))
									{
									break;
									}
								sameSequenceLength++;
								}
							
							if(sameSequenceLength!=pep.length())
								{
								System.err.println("Not Same sequence "+kg.getTranscriptId()+" strand "+kg.getStrand() +" ok-up-to-"+sameSequenceLength);
								System.err.println("P:"+pep.toString()+" "+pep.length());
								System.err.println("Q:"+entrySequence+" "+entrySequence.length());
								if(pep.toString().contains("?"))
									{
									System.err.println("RNA:"+pep.getCodingRNA().toString());
									}
									
								}
							if(sameSequenceLength==0) continue;
							
							nodeList= (NodeList)xpath.evaluate("/uniprot/entry/feature",dom,XPathConstants.NODESET);
							for(int i=0;i< nodeList.getLength();i++)
								{
								final Element feat = (Element)nodeList.item(i);
								final String featureType = feat.getAttribute("type");
								
								if(StringUtils.isBlank(featureType) || feat.equals("chain")) continue;
								if(debug) {
									LOG.debug("featureType:"+featureType);
									}
								final  Element location =(Element)xpath.evaluate("location",feat,XPathConstants.NODE);
								if(location==null) {
									if(debug) {
										LOG.debug("no location ");
										}
									continue;
									}
								int pepStart,pepEnd;
								final Element position = (Element)xpath.evaluate("position",location,XPathConstants.NODE);
								if(position!=null && !StringUtils.isBlank(position.getAttribute("position")))
									{
									pepEnd= Integer.parseInt(position.getAttribute("position"));
									pepStart=pepEnd-1;
									}
								else
									{
									final Element begin = (Element)xpath.evaluate("begin",location,XPathConstants.NODE);
									final Element end = (Element)xpath.evaluate("end",location,XPathConstants.NODE);
									if(begin!=null && end!=null && 
											begin.hasAttribute("position") &&
											end.hasAttribute("position")
											)
										{
										pepStart= Integer.parseInt(begin.getAttribute("position"))-1;
										pepEnd= Integer.parseInt(end.getAttribute("position"))-1;
										}
									else
										{
										if(debug) {
											LOG.debug("cannot find begin/end");
											}
										continue;
										}
									}
								
								if(pepEnd>=sameSequenceLength) {
									if(debug) {
										LOG.debug("pepEnd ("+pepEnd+")>=sameSequenceLength("+sameSequenceLength+")");
										}
									continue;
									}
								
								if(debug) {
									LOG.debug(featureType+" pepStart:"+pepStart+" pepEnd:"+pepEnd+" "+kg.getTranscriptId());
									}
								
								final List<Integer> genomicPos=new ArrayList<Integer>();
								while(pepStart< pepEnd)
									{
									if(pepStart>=pep.length())
										{
										if(debug) {
											System.err.println("pepStart > pep.length()");
											System.err.println("P:"+pep.toString()+" "+pep.length()+" "+kg.getStrand());
											System.err.println("Q:"+entrySequence+" "+ entrySequence.length());
											}
										}
									final int codon[]=pep.convertToGenomic0Coordinates(pepStart);
									pepStart++;
									for(int codon_idx=0;codon_idx<codon.length ;++codon_idx) {
										final int gP = codon[codon_idx];
										if(gP==-1)
											{
											LOG.error("error in genomic for ("+pepStart+"/"+pepEnd+"):");
											System.exit(-1);
											}
										if(debug) {
											System.err.println(kg.getTranscriptId()+" pepStart="+pepStart+" codon["+codon_idx+"] genomic="+gP);
											}
										genomicPos.add(gP);
										}
									}
								Collections.sort(genomicPos);
								final UBed ubed=new UBed();
								ubed.contig = genomic.getContig();
								/** group consecutive coordinates in 'Range' */
								int k0=0;
								while(k0< genomicPos.size())
									{
									Range range=new Range(genomicPos.get(k0), genomicPos.get(k0)+1);
									k0++;
									while(
										k0 < genomicPos.size() && 
										range.end == genomicPos.get(k0)
										)
										{
										range=new Range(range.start,range.end+1);
										++k0;
										}
									ubed.positions.add(range);
									}
								ubed.strand=(byte)(kg.isPositiveStrand()?'+':'-');
								ubed.name=	(String)xpath.evaluate("/uniprot/entry/name",dom,XPathConstants.STRING);
								ubed.featureType = featureType;
								ubed.featureDesc = feat.getAttribute("description");
								ubed.print(this.formatOut,pw);
								}
							}					
							
						}
					pw.flush();
					}//end PrintWriter
				rx.close();
				} // end file reader
			LOG.info("End scan uniprot");
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		return 0;
		}
	
	public static void main(final String[] args)
		{
		new MapUniProtFeatures().instanceMainWithExit(args);
		}
	}
