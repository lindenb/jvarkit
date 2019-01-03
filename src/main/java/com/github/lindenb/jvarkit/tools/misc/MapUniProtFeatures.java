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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLResolver;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;

import org.uniprot.Entry;
import org.uniprot.FeatureType;
import org.uniprot.LocationType;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;

/**

BEGIN_DOC

##Example

```bash
$ java  -jar dist/mapuniprot.jar \
	-R /path/to/human_g1k_v37.fasta \
	-u /path/uri/uniprot.org/uniprot_sprot.xml.gz  \
	-k <(curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk -F '        ' '{if($2 ~ ".*_.*") next; OFS="       "; gsub(/chr/,"",$2);print;}'   ) |\
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
	description="map uniprot features on reference genome",
	keywords={"uniprot","bed","fasta","reference","xjc","xml"}
	)
public class MapUniProtFeatures extends Launcher
	{
	private static final String UNIPROT_NS="http://uniprot.org/uniprot";
	@SuppressWarnings("unused")
	private org.uniprot.ObjectFactory forceJavacCompiler=null;
	private static Logger LOG=Logger.build(MapUniProtFeatures.class).make();
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private Map<String,List<KnownGene>> prot2genes=new HashMap<String,List<KnownGene>>();
    private static class Range
    	{
    	Range(int start,int end)
    		{
    		this.start=start;
    		this.end=end;
    		}
    	int start;
    	int end;
    	}
	private  class UBed
		{
		int tid;
		List<Range> positions=new ArrayList<Range>();
		byte strand;
		String name;
		
		String getChrom()
			{
			return indexedFastaSequenceFile.getSequenceDictionary().getSequence(tid).getSequenceName();
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
		public String toString()
			{
			StringBuilder b=new StringBuilder();
			b.append(getChrom());
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
			return b.toString();
			}
		void print(PrintWriter pw)
			{
			pw.println(this.toString());
			}
		}
	@Parameter(names="-k",description=KnownGene.OPT_KNOWNGENE_DESC)
	private String kgUri=KnownGene.getDefaultUri();
	
	@Parameter(names="-u",description="Uniprot.xml.gz URL/File.")
	private String UNIPROT="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz";
	

	@Parameter(names="-R",description=Launcher.INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File REF=null;
	@Parameter(names="-o",description=OPT_OUPUT_FILE_OR_STDOUT)
	private File OUT=null;

	
	
	
	@Override
	public int doWork(List<String> args)
		{


		
		
		PrintWriter pw=null;
		try
			{
			JAXBContext jc = JAXBContext.newInstance("org.uniprot");
			Unmarshaller uniprotUnmarshaller=jc.createUnmarshaller();
			Marshaller uniprotMarshaller=jc.createMarshaller();

			
			pw=super.openFileOrStdoutAsPrintWriter(OUT);
			LOG.info("read "+REF);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
			if(this.indexedFastaSequenceFile.getSequenceDictionary()==null)
				{
				LOG.error("fasta sequence is not indexed");
				return -1;
				}
			LOG.info("readubf "+kgUri);
			String line;
			Pattern tab=Pattern.compile("[\t]");
			BufferedReader r=IOUtils.openURIForBufferedReading(kgUri);
			while((line=r.readLine())!=null)
				{
				String tokens[]=tab.split(line);
				
				KnownGene kg=new KnownGene();
				kg.setName(tokens[0]);
				kg.setChrom(tokens[1]);
				kg.setStrand(tokens[2].charAt(0));
				kg.setTxStart(Integer.parseInt(tokens[3]));
				kg.setTxEnd(Integer.parseInt(tokens[4]));
				kg.setCdsStart(Integer.parseInt(tokens[5]));
				kg.setCdsEnd(Integer.parseInt(tokens[6]));
				kg.setExonBounds(Integer.parseInt(tokens[7]), tokens[8], tokens[9]);
				List<KnownGene> L=prot2genes.get(tokens[10]);
				if(L==null)
					{
					L=new ArrayList<KnownGene>();
					prot2genes.put(tokens[10], L);
					}
				
				if(indexedFastaSequenceFile.getSequenceDictionary().getSequence(kg.getContig())==null)
					{
					LOG.info("ignoring "+line);
					continue;
					}
				
				L.add(kg);
				
				}
			r.close();
			
			if(OUT!=null) 
				{
				LOG.info("opening "+OUT);
				pw=new PrintWriter(OUT);
				}
			
			LOG.info("knownGenes: "+this.prot2genes.size());
			
			XMLInputFactory xmlInputFactory=XMLInputFactory.newFactory();
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
			
			
			
			LOG.info("Scanning "+UNIPROT);
			
			Reader fr=IOUtils.openURIForBufferedReading(UNIPROT);
			XMLEventReader rx=xmlInputFactory.createXMLEventReader(fr);
		
			final QName uEntry=new QName(UNIPROT_NS,"entry");
			GenomicSequence genomic=null;
			while(rx.hasNext())
				{
				
				XMLEvent evt=rx.peek();
				if(!(evt.isStartElement() && evt.asStartElement().getName().equals(uEntry)))
					{
					rx.next();
					continue;
					}
				JAXBElement<Entry> jaxbElement=uniprotUnmarshaller.unmarshal(rx, Entry.class);
				Entry entry= jaxbElement.getValue();
				
				

				
				if(entry.getFeature().isEmpty()) continue;
				List<KnownGene> genes=null;
				
				for(String acn:entry.getAccession())
					{
					genes=prot2genes.get(acn);
					if(genes!=null) break;
					}
				if(genes==null) continue;
				
				for(KnownGene kg:genes)
					{
					if(genomic==null ||  !genomic.getChrom().equals(kg.getChromosome()))
						{
						genomic=new GenomicSequence(this.indexedFastaSequenceFile, kg.getChromosome());
						}
					
					
					KnownGene.Peptide pep=kg.getCodingRNA(genomic).getPeptide();
					
					
					int sameSequenceLength=0;
					while(  sameSequenceLength < entry.getSequence().getValue().length() &&
							sameSequenceLength < pep.length() 
							)
						{
						if(Character.toUpperCase(entry.getSequence().getValue().charAt(sameSequenceLength))!=Character.toUpperCase(entry.getSequence().getValue().charAt(sameSequenceLength)))
							{
							break;
							}
						sameSequenceLength++;
						}
					
					if(sameSequenceLength!=pep.length())
						{
						System.err.println("Not Same sequence "+kg.getName()+" strand "+kg.getStrand() +" ok-up-to-"+sameSequenceLength);
						System.err.println("P:"+pep.toString()+" "+pep.length());
						System.err.println("Q:"+entry.getSequence().getValue()+" "+entry.getSequence().getLength());
						if(pep.toString().contains("?"))
							{
							System.err.println("RNA:"+pep.getCodingRNA().toString());
							
							}
							
						}
					if(sameSequenceLength==0) continue;
					
					for(FeatureType feat:entry.getFeature())
						{
						if(feat.getType()==null || feat.getType().isEmpty()) continue;
						LocationType locType=feat.getLocation();
						if(locType==null) continue;
						int pepStart,pepEnd;
						if(locType.getPosition()!=null && locType.getPosition().getPosition()!=null)
							{
							pepEnd=locType.getPosition().getPosition().intValue();
							pepStart=pepEnd-1;
							}
						else if(locType.getBegin()!=null &&
								locType.getEnd()!=null &&
								locType.getBegin().getPosition()!=null &&
								locType.getEnd().getPosition()!=null )
							{
							pepStart=locType.getBegin().getPosition().intValue()-1;
							pepEnd=locType.getEnd().getPosition().intValue();
							}
						else
							{
							continue;
							}
						if(pepEnd>=sameSequenceLength) continue;
						
						List<Integer> genomicPos=new ArrayList<Integer>();
						while(pepStart< pepEnd)
							{
							if(pepStart>=pep.length())
								{
								System.err.println();
								System.err.println("P:"+pep.toString()+" "+pep.length()+" "+kg.getStrand());
								System.err.println("Q:"+entry.getSequence().getValue()+" "+entry.getSequence().getLength());
								uniprotMarshaller.marshal(new JAXBElement<FeatureType>(new QName(UNIPROT_NS, "feature"), FeatureType.class,feat),System.err);

								}
							int codon[]=pep.convertToGenomicCoordinates(pepStart);
							pepStart++;
							for(int gP:codon)
								{
								if(gP==-1)
									{
									uniprotMarshaller.marshal(new JAXBElement<FeatureType>(new QName(UNIPROT_NS, "feature"), FeatureType.class,feat),System.err);
									LOG.error("error in genomoc for ("+pepStart+"/"+pepEnd+"):"+entry.getName());
									System.exit(-1);
									}
								genomicPos.add(gP);
								}
							}
						Collections.sort(genomicPos);
						UBed ubed=new UBed();
						ubed.tid=this.indexedFastaSequenceFile.getSequenceDictionary().getSequenceIndex(genomic.getChrom());
						int k0=0;
						while(k0< genomicPos.size())
							{
							Range range=new Range(genomicPos.get(k0), genomicPos.get(k0)+1);
							k0++;
							while(
								k0< genomicPos.size() && 
								range.end==genomicPos.get(k0)
								)
								{
								range=new Range(range.start,range.end+1);
								++k0;
								}
							ubed.positions.add(range);
							}
						ubed.strand=(byte)(kg.isPositiveStrand()?'+':'-');
						ubed.name=feat.getType();
						ubed.print(pw);
						}
					}					
					
				}
			rx.close();
			fr.close();
			LOG.info("End scan uniprot");
/*
			mappedFeatures.doneAdding();
			
			
			CloseableIterator<UBed> iter=mappedFeatures.iterator();
			while(iter.hasNext())
				{
				UBed ubed=iter.next();
				ubed.print();
				}
			mappedFeatures.cleanup();
			*/
			}
		catch(Exception err)
			{
			err.printStackTrace();
			if(OUT!=null) OUT.deleteOnExit();
			return -1;
			}
		finally
			{
			pw.flush();
			pw.close();
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		return 0;
		}
	
	public static void main(String[] args)
		{
		new MapUniProtFeatures().instanceMainWithExit(args);
		}
	}
