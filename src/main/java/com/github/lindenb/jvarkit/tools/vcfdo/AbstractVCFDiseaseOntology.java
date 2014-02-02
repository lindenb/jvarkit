package com.github.lindenb.jvarkit.tools.vcfdo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.Reader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broadinstitute.variant.variantcontext.VariantContext;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;



import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.Function;
import com.github.lindenb.jvarkit.util.biomart.BiomartQuery;
import com.github.lindenb.jvarkit.util.cli.GetOpt;
import com.github.lindenb.jvarkit.util.doid.DiseaseOntoglogyTree;
import com.github.lindenb.jvarkit.util.picard.IntervalTreeMapFactory;
import com.github.lindenb.jvarkit.util.picard.SamSequenceRecordTreeMap;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;

public abstract class AbstractVCFDiseaseOntology
	extends AbstractVCFFilter2
	{
	/* Disease Ontology OWL file/URI. */
	protected String DOI_INPUT="http://www.berkeleybop.org/ontologies/doid.owl";
	/* Disease Ontology Annotations.",optional=true */
	protected String DOI_ANN="http://dga.nubic.northwestern.edu/ajax/Download.ajax.php?exportType=ids";
	/*  used to reduce the number of mapped genes */
	protected File REF=null;
	
	
	protected DiseaseOntoglogyTree diseaseOntoglogyTree;
	protected Map<Integer,Set<DiseaseOntoglogyTree.Term>> gene2doid=new HashMap<Integer,Set<DiseaseOntoglogyTree.Term>>();
	protected Map<String,Set<DiseaseOntoglogyTree.Term>> ensemblProtein2doid=new HashMap<String,Set<DiseaseOntoglogyTree.Term>>();
	private SamSequenceRecordTreeMap<Integer> ncbiGeneMap=null;
	@SuppressWarnings("unused")
	private SamSequenceRecordTreeMap<String> ensemblProteinMap;
	
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -D (url) Disease Ontology OWL file/URI. default:"+DOI_INPUT);
		out.println(" -A (url) Disease Annotations. default:"+DOI_ANN);
		out.println(" -R (fasta) Indexed  Genome Reference.");
		super.printOptions(out);
		}
	
	@Override
	protected String getGetOptDefault()
		{
		return super.getGetOptDefault()+"D:A:R:";
		}
	@Override
	protected GetOptStatus handleOtherOptions(int c, GetOpt opt, String[] args)
		{
		switch(c)
			{
			case 'D': DOI_INPUT=opt.getOptArg(); return GetOptStatus.OK;
			case 'A': DOI_ANN=opt.getOptArg(); return GetOptStatus.OK;
			case 'R': REF=new File(opt.getOptArg());return GetOptStatus.OK;
			default:return super.handleOtherOptions(c, opt, args);
			}
		
		}
	
	protected void readDiseaseOntoglogyTree() throws IOException
		{
		this.info("read DOI "+DOI_INPUT);
		try
			{
			diseaseOntoglogyTree=DiseaseOntoglogyTree.parse(DOI_INPUT);
			this.info("GO size:"+diseaseOntoglogyTree.size());
			}
		catch(XMLStreamException err)
			{
			throw new IOException(err);
			}
		}
	protected void readJensenLabAnnotations() throws IOException
		{
		this.info("read DOID-Annotation from "+DOI_ANN);
		BufferedReader in=IOUtils.openURIForBufferedReading(DOI_ANN);
		String line;
		while((line=in.readLine())!=null)
			{
			if(!line.startsWith("ENSP")) continue;
			}
		in.close();
		}
	
	protected void readDiseaseOntoglogyAnnotations() throws IOException
		{
		try
			{
			this.info("read DOID-Annotation from "+DOI_ANN);
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
	
	protected void loadEnsemblProtein() throws IOException
		{
	
		BiomartQuery q=new BiomartQuery();
		q.setDataSetName("hsapiens_gene_ensembl");
		q.setAttributes(
				"chromosome_name",
				"start_position",
				"end_position",
				"ensembl_peptide_id"
				);
		q.setUniqRows(true);
		this.info("sending "+q);
		
		IntervalTreeMapFactory<String> itf=new IntervalTreeMapFactory<String>();
		if(REF!=null)
			{
			itf.setSamSequenceDictionary(new IndexedFastaSequenceFile(REF).getSequenceDictionary());
			}
		itf.setValueFunction(new Function<String[], String>()
			{
			@Override
			public String apply(String[] param)
				{
				if(param.length<4 || param[3].isEmpty()) return null;
				String ensp=param[3];
				
				if(!ensemblProtein2doid.containsKey(ensp)) return null;
				return ensp;
				}
			});
		this.info("invoking biomart "+q);
		LineReader r=q.execute();
		
		this.ensemblProteinMap=itf.createIntervalMap(r);
		r.close();
		}
	
	
	protected void loadEntrezGenes(SAMSequenceDictionary dict) throws IOException
		{
		this.ncbiGeneMap=new SamSequenceRecordTreeMap<Integer>(dict);
		BiomartQuery q=new BiomartQuery();
		q.setDataSetName("hsapiens_gene_ensembl");
		q.setAttributes(
				"chromosome_name",
				"start_position",
				"end_position",
				"entrezgene"
				);
		q.setUniqRows(true);
		this.info("sending "+q);
		
		
		this.info("invoking biomart "+q);
		LineReader r=q.execute();
		LineIterator lr=new LineIteratorImpl(r);
		Pattern pattern=Pattern.compile("[\t]");
		while(lr.hasNext())
			{
			String line=lr.next();
			String param[]=pattern.split(line);
			if(param.length<4 || param[3].isEmpty()) continue;
			SAMSequenceRecord rec=null;
			if((rec=dict.getSequence(param[0]))==null) continue;
			try
				{
				int start=Integer.parseInt(param[0]);
				int end=Integer.parseInt(param[1]);
				int ncbiGene=Integer.parseInt(param[3]);
				if(!gene2doid.containsKey(ncbiGene)) continue;
				ncbiGeneMap.put(rec.getSequenceIndex(),start,end,ncbiGene);
				}
			catch (Exception e)
				{
				warning(e);
				continue;
				}
			
			
			}
		r.close();
		}

	protected Set<Integer> getGeneIds(VariantContext ctx)
		{
		return new HashSet<Integer>(this.ncbiGeneMap.getOverlapping(
				ctx.getChr(),
				ctx.getStart(),
				ctx.getEnd()
				));
		}
	
	}
