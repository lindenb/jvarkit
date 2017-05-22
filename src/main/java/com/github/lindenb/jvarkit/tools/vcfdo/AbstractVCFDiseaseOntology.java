
/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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

*/package com.github.lindenb.jvarkit.tools.vcfdo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
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

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.biomart.BiomartQuery;
import com.github.lindenb.jvarkit.util.doid.DiseaseOntoglogyTree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

public abstract class AbstractVCFDiseaseOntology
	extends Launcher
	{
	private static final Logger LOG= Logger.build(AbstractVCFDiseaseOntology.class).make();
	
	/* Disease Ontology OWL file/URI. */
	@Parameter(names="-D",description="Disease Ontology OWL file/URI")
	protected String DOI_INPUT="http://www.berkeleybop.org/ontologies/doid.owl";
	/* Disease Ontology Annotations.",optional=true */
	@Parameter(names="-A",description="Disease Annotations")
	protected String DOI_ANN="http://dga.nubic.northwestern.edu/ajax/Download.ajax.php?exportType=ids";
	/*  used to reduce the number of mapped genes */
	@Parameter(names="-R",description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	protected File REF=null;
	
	
	protected DiseaseOntoglogyTree diseaseOntoglogyTree;
	protected Map<Integer,Set<DiseaseOntoglogyTree.Term>> gene2doid=new HashMap<Integer,Set<DiseaseOntoglogyTree.Term>>();
	protected Map<String,Set<DiseaseOntoglogyTree.Term>> ensemblProtein2doid=new HashMap<String,Set<DiseaseOntoglogyTree.Term>>();
	private IntervalTreeMap<Integer> ncbiGeneMap=null;
	private IntervalTreeMap<String> ensemblProteinMap;
	
	
	
	
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
	protected void readJensenLabAnnotations() throws IOException
		{
		LOG.info("read DOID-Annotation from "+DOI_ANN);
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
		LOG.info("sending "+q);
		LOG.info("invoking biomart "+q);
		LineReader r=q.execute();
		LineIteratorImpl iter = new LineIteratorImpl(r);
		while(iter.hasNext()) {
			String line=iter.next();
			String param[]=line.split("[\t]");
			if(param.length<4 || param[3].isEmpty()) continue;
			String ensp=param[3];
			
			if(!ensemblProtein2doid.containsKey(ensp)) continue;
			ensemblProteinMap.put(new Interval(param[0],Integer.parseInt(param[1]),Integer.parseInt(param[2])),ensp);
		}
		
		iter.close();
		r.close();
		}
	
	
	protected void loadEntrezGenes(SAMSequenceDictionary dict) throws IOException
		{
		this.ncbiGeneMap=new IntervalTreeMap<Integer>();
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
		
		
		LOG.info("invoking biomart "+q);
		LineReader r=q.execute();
		@SuppressWarnings("resource")
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
				ncbiGeneMap.put(new Interval(rec.getSequenceName(),start,end),ncbiGene);
				}
			catch (Exception e)
				{
				LOG.warning(e);
				continue;
				}
			
			
			}
		r.close();
		}

	protected Set<Integer> getGeneIds(VariantContext ctx)
		{
		return new HashSet<Integer>(this.ncbiGeneMap.getOverlapping(new Interval(
				ctx.getContig(),
				ctx.getStart(),
				ctx.getEnd()
				)));
		}
	
	}
