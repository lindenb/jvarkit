package com.github.lindenb.jvarkit.tools.vcfdo;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;

import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.picard.util.Log;


import com.github.lindenb.jvarkit.lang.Function;
import com.github.lindenb.jvarkit.util.AbstractVCFFilter;
import com.github.lindenb.jvarkit.util.biomart.BiomartQuery;
import com.github.lindenb.jvarkit.util.doid.DiseaseOntoglogyTree;
import com.github.lindenb.jvarkit.util.picard.IOUtils;
import com.github.lindenb.jvarkit.util.picard.IntervalTreeMapFactory;

public abstract class AbstractVCFDiseaseOntology extends AbstractVCFFilter
	{
	@Option(shortName="DOI", doc="Disease Ontology OWL file/URI.",optional=true)
	public String DOI_INPUT="http://www.berkeleybop.org/ontologies/doid.owl";
	@Option(shortName="DGA", doc="Disease Ontology Annotations.",optional=true)
	public String DOI_ANN="http://dga.nubic.northwestern.edu/ajax/Download.ajax.php?exportType=ids";
	@Option(shortName=StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="Reference (used to reduce the number of mapped genes)",optional=true)
	public File REF=null;

	protected static final Log LOG=Log.getInstance(AbstractVCFDiseaseOntology.class);

	
	protected DiseaseOntoglogyTree diseaseOntoglogyTree;
	protected Map<Integer,Set<DiseaseOntoglogyTree.Term>> gene2doid=new HashMap<Integer,Set<DiseaseOntoglogyTree.Term>>();
	private IntervalTreeMap<Integer> ncbiGeneMap;

	
	protected void readDiseaseOntoglogyTree() throws IOException
		{
		LOG.info("read DOI "+DOI_INPUT);
		try
			{
			diseaseOntoglogyTree=DiseaseOntoglogyTree.parse(DOI_INPUT);
			LOG.info("GO size:"+diseaseOntoglogyTree.size());
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	
	protected void readDiseaseOntoglogyAnnotations() throws IOException
		{
		try
			{
			LOG.info("read DOID-Annotation from "+DOI_ANN);
		XMLInputFactory xmlInputFactory = XMLInputFactory.newInstance();
		xmlInputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
		xmlInputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
		xmlInputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
		Reader in=IOUtils.openURIForBufferedReading(DOI_ANN);
		XMLEventReader reader= xmlInputFactory.createXMLEventReader(in);	
		final String RDFNS="http://www.w3.org/1999/02/22-rdf-syntax-ns#";
		final String CDNS="http://www.ncbi.nlm.nih.gov/";
		final QName descriptionQName=new QName(RDFNS, "Description");
		while(reader.hasNext())
			{
			XMLEvent evt=reader.nextEvent();
			if(!evt.isStartElement()) continue;
			StartElement E1=evt.asStartElement();
			if(!E1.getName().equals(descriptionQName)) continue;
			Set<Integer> genes=new HashSet<Integer>();
			Set<DiseaseOntoglogyTree.Term> terms=new HashSet<DiseaseOntoglogyTree.Term>();
			while(reader.hasNext())
				{
				evt=reader.nextEvent();
				if(evt.isStartElement())
					{
					StartElement E2=evt.asStartElement();
					if(!CDNS.equals(E2.getName().getNamespaceURI())) continue;
					if(E2.getName().equals("DOID"))
						{
						String doid=reader.getElementText();
						if(doid.startsWith("DOID:")) doid="DOID:"+doid;
						DiseaseOntoglogyTree.Term t=diseaseOntoglogyTree.getTermByAccession(doid);
						if(t!=null)
							{
							terms.add(t);
							}
						}
					else if(E2.getName().equals("GeneId"))
						{
						genes.add(Integer.parseInt(reader.getElementText()));
						}
					}
				else if(evt.isEndElement())
					{
					if(evt.asEndElement().getName().equals(descriptionQName))
						{
						for(Integer geneId:genes)
							{
							Set<DiseaseOntoglogyTree.Term> set=gene2doid.get(geneId);
							if(set==null)
								{
								set=new HashSet<DiseaseOntoglogyTree.Term>();
								gene2doid.put(geneId,set);
								}
							set.addAll(terms);
							}
						
						break;
						}
					}
				}
			
			}
		in.close();
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	protected void loadEntrezGenes() throws IOException
		{
	
		BiomartQuery q=new BiomartQuery();
		q.setDataSetName("hsapiens_gene_ensembl");
		q.setAttributes(
				"chromosome_name",
				"start_position",
				"end_position",
				"entrezgene"
				);
		q.setUniqRows(true);
		LOG.info("sending "+q);
		
		IntervalTreeMapFactory<Integer> itf=new IntervalTreeMapFactory<Integer>();
		if(REF!=null)
			{
			itf.setSamSequenceDictionary(new IndexedFastaSequenceFile(REF).getSequenceDictionary());
			}
		itf.setValueFunction(new Function<String[], Integer>()
			{
			@Override
			public Integer apply(String[] param)
				{
				if(param.length<4 || param[3].isEmpty()) return null;
				Integer ncbiGene=null;
				try
					{
					ncbiGene=new Integer(param[3]);
					}
				catch (Exception e)
					{
					return null;
					}
				if(!gene2doid.containsKey(ncbiGene)) return null;
				return ncbiGene;
				}
			});
		LOG.info("invoking biomart "+q);
		LineReader r=q.execute();
		
		this.ncbiGeneMap=itf.createIntervalMap(r);
		r.close();
		}

	protected Set<Integer> getGeneIds(VariantContext ctx)
		{
		return new HashSet<Integer>(this.ncbiGeneMap.getOverlapping(new Interval(
				ctx.getChr(),
				ctx.getStart(),
				ctx.getEnd()
				)));
		}

	}
