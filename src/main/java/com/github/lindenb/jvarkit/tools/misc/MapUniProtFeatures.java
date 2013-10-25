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

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.illumina.parser.Range;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Log;
import net.sf.picard.util.ProgressLogger;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.ucsc.KnownGene;
import org.uniprot.*;

@SuppressWarnings("unused")
public class MapUniProtFeatures extends AbstractCommandLineProgram
	{
	private static final String UNIPROT_NS="http://uniprot.org/uniprot";
	private org.uniprot.ObjectFactory forceJavacCompiler=null;
	private static Log LOG=Log.getInstance(MapUniProtFeatures.class);
	
    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " map uniprot features on reference genome. ";

	
    @Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME,doc="Reference",optional=false)
    public File REF;
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="output name (default: stdout)",optional=false)
    public File OUT=null;

	
	
    @Option(shortName="KG",doc="KnownGene data URI/File. should look like http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz . Beware chromosome names are formatted the same as your REFERENCE.",optional=false)
	public String kgUri="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz";
	
	private Map<String,List<KnownGene>> prot2genes=new HashMap<String,List<KnownGene>>();
	
    @Option(shortName="UP",doc="Uniprot URL/File",optional=true)
	public String UNIPROT="ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz";

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
    
	@Override
	protected int doWork()
		{
		PrintWriter pw=new PrintWriter(System.out);
		try
			{
			JAXBContext jc = JAXBContext.newInstance("org.uniprot");
			Unmarshaller uniprotUnmarshaller=jc.createUnmarshaller();
			Marshaller uniprotMarshaller=jc.createMarshaller();

			
			LOG.info("read "+REF);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(REF);
			LOG.info("readubf "+kgUri);
			String line;
			Pattern tab=Pattern.compile("[\t]");
			BufferedReader r=IOUtils.openURIForBufferedReading(this.kgUri);
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
				
				if(indexedFastaSequenceFile.getSequenceDictionary().getSequence(kg.getChr())==null)
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
						if(super.VALIDATION_STRINGENCY!=SAMFileReader.ValidationStringency.SILENT )
							{
							System.err.println("Not Same sequence "+kg.getName()+" strand "+kg.getStrand() +" ok-up-to-"+sameSequenceLength);
							System.err.println("P:"+pep.toString()+" "+pep.length());
							System.err.println("Q:"+entry.getSequence().getValue()+" "+entry.getSequence().getLength());
							if(pep.toString().contains("?"))
								{
								System.err.println("RNA:"+pep.getCodingRNA().toString());
								System.exit(-1);
								}
							
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
