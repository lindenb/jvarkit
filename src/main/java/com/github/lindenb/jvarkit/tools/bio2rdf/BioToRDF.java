package com.github.lindenb.jvarkit.tools.bio2rdf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.xml.namespace.QName;
import javax.xml.stream.EventFilter;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import org.apache.jena.rdf.model.Model;
import org.apache.jena.rdf.model.ModelFactory;
import org.apache.jena.rdf.model.Property;
import org.apache.jena.rdf.model.Resource;
import org.apache.jena.rdf.model.ResourceFactory;
import org.apache.jena.rdf.model.Statement;
import org.apache.jena.rdfxml.xmlinput.StAX2Model;
import org.apache.jena.util.iterator.ExtendedIterator;
import org.apache.jena.vocabulary.DC;
import org.apache.jena.vocabulary.OWL;
import org.apache.jena.vocabulary.RDF;
import org.apache.jena.vocabulary.RDFS;
import org.apache.jena.vocabulary.XSD;
import org.xml.sax.SAXException;

import com.beust.jcommander.DynamicParameter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jena.JenaUtils;
import com.github.lindenb.jvarkit.jena.vocabulary.OBOInOwl;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.iterator.LineIterators;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.xml.SimpleEventFilter;

import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

/**
BEGIN_DOC

# motivation

build  RDF file to build a RDF database.
 
 
# usage
 
```
$ java -jar dist/jvarkit.jar bio2rdf | gzip > bio2rdf.rdf.gz
```

## Example queries:

### Example

find all the descendant of `GO:0045823`

```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>


SELECT ?subClass ?label ?desc WHERE { 
        ?subClass rdfs:subClassOf* <http://purl.obolibrary.org/obo/GO_0045823> . 
        ?subClass dc:title ?label . 
        ?subClass dc:description ?desc . 
    }
```

```
$ arq --data=/home/me/bio2rdf.rdf.gz --query query.sparql
```

| subClass                                    | label        | desc                                                                                                                                    |
|--- |-- |-- |
| <http://purl.obolibrary.org/obo/GO_0045823> | "GO:0045823" | "positive regulation of heart contraction"                                                                                              |
| <http://purl.obolibrary.org/obo/GO_0001989> | "GO:0001989" | "positive regulation of the force of heart contraction involved in baroreceptor response to decreased systemic arterial blood pressure" |
| <http://purl.obolibrary.org/obo/GO_0010460> | "GO:0010460" | "positive regulation of heart rate"                                                                                                     |
| <http://purl.obolibrary.org/obo/GO_0001996> | "GO:0001996" | "positive regulation of heart rate by epinephrine-norepinephrine"                                                                       |
| <http://purl.obolibrary.org/obo/GO_0086024> | "GO:0086024" | "adenylate cyclase-activating adrenergic receptor signaling pathway involved in positive regulation of heart rate"                      |
| <http://purl.obolibrary.org/obo/GO_0003065> | "GO:0003065" | "positive regulation of heart rate by epinephrine"                                                                                      |
| <http://purl.obolibrary.org/obo/GO_0003112> | "GO:0003112" | "positive regulation of heart rate by neuronal epinephrine"                                                                             |
| <http://purl.obolibrary.org/obo/GO_0003111> | "GO:0003111" | "positive regulation of heart rate by circulating epinephrine"                                                                          |
| <http://purl.obolibrary.org/obo/GO_0003066> | "GO:0003066" | "positive regulation of heart rate by norepinephrine"                                                                                   |
| <http://purl.obolibrary.org/obo/GO_0003114> | "GO:0003114" | "positive regulation of heart rate by circulating norepinephrine"                                                                       |
| <http://purl.obolibrary.org/obo/GO_0003113> | "GO:0003113" | "positive regulation of heart rate by neuronal norepinephrine"                                                                          |
| <http://purl.obolibrary.org/obo/GO_0001988> | "GO:0001988" | "positive regulation of heart rate involved in baroreceptor response to decreased systemic arterial blood pressure"                     |
| <http://purl.obolibrary.org/obo/GO_0060452> | "GO:0060452" | "positive regulation of cardiac muscle contraction"                                                                                     |
| <http://purl.obolibrary.org/obo/GO_0003099> | "GO:0003099" | "positive regulation of the force of heart contraction by chemical signal"                                                              |
| <http://purl.obolibrary.org/obo/GO_0003061> | "GO:0003061" | "positive regulation of the force of heart contraction by norepinephrine"                                                               |
| <http://purl.obolibrary.org/obo/GO_0003110> | "GO:0003110" | "positive regulation of the force of heart contraction by neuronal norepinephrine"                                                      |
| <http://purl.obolibrary.org/obo/GO_0003109> | "GO:0003109" | "positive regulation of the force of heart contraction by circulating norepinephrine"                                                   |
| <http://purl.obolibrary.org/obo/GO_0003059> | "GO:0003059" | "positive regulation of the force of heart contraction by epinephrine"                                                                  |
| <http://purl.obolibrary.org/obo/GO_0003087> | "GO:0003087" | "positive regulation of the force of heart contraction by neuronal epinephrine"                                                         |
| <http://purl.obolibrary.org/obo/GO_0003088> | "GO:0003088" | "positive regulation of the force of heart contraction by circulating epinephrine"                                                      |
| <http://purl.obolibrary.org/obo/GO_0001997> | "GO:0001997" | "positive regulation of the force of heart contraction by epinephrine-norepinephrine"                                                   |
| <http://purl.obolibrary.org/obo/GO_0003090> | "GO:0003090" | "positive regulation of the force of heart contraction by neuronal epinephrine-norepinephrine"                                          |
| <http://purl.obolibrary.org/obo/GO_0003089> | "GO:0003089" | "positive regulation of the force of heart contraction by circulating epinephrine-norepinephrine"                                       |

### find all the phenotypes containing the word 'arrhythmia'

```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX bio: <https://umr1087.univ-nantes.fr/bio2rdf/>

SELECT  ?label ?desc WHERE { 
        ?disease a bio:Phenotype .
        ?disease dc:title ?label . 
        ?disease dc:description ?desc . 
	FILTER (CONTAINS(LCASE(?desc),"arrhythmia")) .
    }
```

```
------------------------------------------------
| label        | desc                          |
================================================
| "HP:0002521" | "Hypsarrhythmia"              |
| "HP:0011215" | "Hemihypsarrhythmia"          |
| "HP:0004308" | "Ventricular arrhythmia"      |
| "HP:0001692" | "Atrial arrhythmia"           |
| "HP:0005115" | "Supraventricular arrhythmia" |
| "HP:0011675" | "Arrhythmia"                  |
------------------------------------------------
```


END_DOC
 
**/
@Program(
		name="bio2rdf",
		description="Build a RDF database for human from misc sources",
		keywords={"rdf","ontology","sparql"},
		creationDate="20220427",
		modificationDate="20220427",
		jvarkit_amalgamion =  true
		)
public class BioToRDF extends Launcher {
	private static final Logger LOG = Logger.build(BioToRDF.class).make();
	private final static String PREFIX="bio";
	private final static String NS = "https://umr1087.univ-nantes.fr/bio2rdf/";
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"--genes"},description="Limit to those genes names , separated with comma (for debugging)")
	private String limitGenesStr="";

	
	
	@SuppressWarnings("serial")
	@DynamicParameter(names={"-D"},description="parameters. -Dkey1=value1  -Dkey2=value2 ...")
	private Map<String,String> resourceMap = new HashMap<String,String>() {{{
		put("NCBI_GENE_INFO", "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz");
		put("NCBI_GENE_GO","https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz");
		put("GENCODE_RELEASE", "43");
		put("GFF3_GRCH38", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.annotation.gff3.gz");
		put("GFF3_GRCH37", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/GRCh37_mapping/gencode.v{GENCODE_RELEASE}lift37.annotation.gff3.gz");
		//put("HUMAN_DO_OWL", "https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/main/src/ontology/HumanDO.owl");
		put("HUMAN_HPO_OWL", "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/hp.owl");
		put("HPO_PHENOTYPE_TO_GENE", "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/phenotype_to_genes.txt");
		put("GO_OWL","http://purl.obolibrary.org/obo/go.owl");
		put("BIOGRID_XML_25","https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.221/BIOGRID-ALL-4.4.221.psi25.zip");
		//put("MONDO_OWL","https://github.com/monarch-initiative/mondo/releases/download/v2023-04-04/mondo.owl");
		}}};


	private XMLStreamWriter writer = null;
	private final Map<String, NcbiGeneInfo> geneid2gene = new HashMap<>();
	private final Map<String, NcbiGeneInfo> symbol2gene = new HashMap<>();
	private final Map<String, NcbiGeneInfo> hgnc2gene = new HashMap<>();
	
	private static class NcbiGeneInfo {
		String geneid = null;
		String symbol = null;
		String hgnc = null;
		public String getURI() {
			if(hgnc!=null && hgnc.startsWith("HGNC:")) 
				{
				final String hgnc_id = this.hgnc.substring(5);
				return "http://identifiers.org/hgnc/"+ hgnc_id;
				}
			return "https://www.ncbi.nlm.nih.gov/gene/"+ this.geneid;
			}
		@Override
		public int hashCode() {
			return geneid.hashCode();
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof NcbiGeneInfo)) return false;
			return this.geneid.equals(NcbiGeneInfo.class.cast(obj).geneid);
			}
		}
	
	private void parseHPOA(String uri) throws IOException,XMLStreamException {
		LOG.info("parsing "+uri);
		String line;
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			while((line=br.readLine())!=null) {
				if(line.startsWith("#")) continue;
				break;
				}
			if(line==null) throw new IOException("Cannot read header line of "+uri);
			final FileHeader header = new FileHeader(CharSplitter.TAB.split(line));
			while((line=br.readLine())!=null) {
				final Map<String,String> rec = header.toMap(CharSplitter.TAB.split(line));
				writer.writeStartElement(PREFIX,"Phenotype",NS);

				writer.writeEndElement();
				}
			}
		}
	
	private void parseNcbiGeneInfo(String uri) throws IOException,XMLStreamException {
		LOG.info("parsing "+uri);
		
		final Set<String> limit_gene_names = Arrays.stream(limitGenesStr.split("[, \t]")).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toSet());
				
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			String line = br.readLine();
			if(line==null || !line.startsWith("#")) throw new IOException("Cannot read first line or "+uri);
			line=line.substring(1);//remove '#'
			final FileHeader header = new FileHeader(CharSplitter.TAB.split(line));
			while((line=br.readLine())!=null) {
				final Map<String,String> rec = header.toMap(CharSplitter.TAB.split(line));
				if(!rec.get("tax_id").equals("9606")) continue;
				//limit gene name for debugging
				if(!limit_gene_names.isEmpty() && !limit_gene_names.contains(rec.get("Symbol"))) {
					continue;
					}
				
				final NcbiGeneInfo info = new NcbiGeneInfo();
				info.geneid = rec.get("GeneID");
				info.symbol = rec.get("Symbol");
				this.geneid2gene.put(info.geneid, info);
				this.symbol2gene.put(info.symbol, info);
				
				writer.writeStartElement(PREFIX,"Gene",NS);
				writer.writeAttribute("rdf", RDF.getURI(), "about", info.getURI() );
				
				writer.writeStartElement("dc","title",DC.getURI());
				writer.writeCharacters(info.symbol);
				writer.writeEndElement();
				
				writer.writeStartElement(PREFIX,"symbol",NS);
				writer.writeCharacters(info.symbol);
				writer.writeEndElement();
				writer.writeStartElement(PREFIX,"geneid",NS);
				writer.writeCharacters(info.geneid);
				writer.writeEndElement();

				
				String s=rec.get("description");
				if(!StringUtils.isBlank(s)) {
					writer.writeStartElement("dc","description",DC.getURI());
					writer.writeCharacters(s);
					writer.writeEndElement();
					}
				
				s=rec.get("type_of_gene");
				if(!StringUtils.isBlank(s) && !s.equals("-")) {
					writer.writeStartElement(PREFIX,"gene_type",NS);
					writer.writeCharacters(s.replace('-', '_'));
					writer.writeEndElement();
					}
				
 				for(String syn: CharSplitter.PIPE.split(rec.get("Synonyms"))) {
 					if(StringUtils.isBlank(syn)) continue;
 					if(syn.equals("-")) continue;
 					writer.writeStartElement(PREFIX,"synonym",NS);
					writer.writeCharacters(syn);
					writer.writeEndElement();
 					}
 				for(String xref: CharSplitter.PIPE.split(rec.get("dbXrefs"))) {
 					if(StringUtils.isBlank(xref)) continue;
 					final int colon = xref.indexOf(':');
 					if(colon==-1) continue;
 					final String key = xref.substring(0,colon);
 					if(key.equals("Ensembl")) {
 						writer.writeStartElement(PREFIX,"ensembl_geneid",NS);
 						writer.writeCharacters(xref.substring(colon+1));
 						writer.writeEndElement();
 						}
 					else if(key.equals("HGNC")) {
 						final String hgnc_id= xref.substring(colon+1);
 						writer.writeStartElement(PREFIX,"hgnc_id",NS);
 						writer.writeCharacters(hgnc_id);
 						writer.writeEndElement();
 						info.hgnc = hgnc_id;
 						hgnc2gene.put(hgnc_id, info);
 						}
 					}
				writer.writeEndElement();
				writer.writeCharacters("\n");
				}
			}
		}
	
	private void parseHPO2gene(String uri) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uri)) return;
		LOG.info("parsing "+uri);
		
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			String line = br.readLine();
			if(line==null) throw new IOException("Cannot read first line or "+uri);
			final FileHeader header = new FileHeader(CharSplitter.TAB.split(line));
			String prev_hpo_id=null;
			final Set<String> ncbi_gene_ids= new HashSet<>();
			for(;;) {
				line=br.readLine();
				final Map<String,String> rec = (line==null?null:header.toMap(CharSplitter.TAB.split(line)));
				if(rec==null || (prev_hpo_id!=null && !prev_hpo_id.equals(rec.get("hpo_id")))) {
					if(prev_hpo_id!=null && 
						!ncbi_gene_ids.isEmpty() &&
						ncbi_gene_ids.stream().anyMatch(ID->geneid2gene.containsKey(ID))) {
						writer.writeStartElement("rdf", "Description", RDF.getURI() );
						writer.writeAttribute("rdf", RDF.getURI(), "about","http://purl.obolibrary.org/obo/" + prev_hpo_id.replace(':', '_'));
						for(String gene_id : ncbi_gene_ids) {
							final NcbiGeneInfo gene = this.geneid2gene.get(gene_id);
							if(gene==null) continue;
							writer.writeEmptyElement(PREFIX, "has_gene", NS);
							writer.writeAttribute("rdf",RDF.getURI(),"resource",gene.getURI());
							}
						writer.writeEndElement();
						writer.writeCharacters("\n");
						}
					
					if(rec==null) break;
					ncbi_gene_ids.clear();
					}
				prev_hpo_id = rec.get("hpo_id");
				ncbi_gene_ids.add(rec.get("ncbi_gene_id"));
				}
			}
		}
	
	private void parseNcbiGeneGO(String uri) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uri)) return;
		LOG.info("parsing "+uri);
		final List<Map.Entry<NcbiGeneInfo,String>> gene2go = new ArrayList<>();
		
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			String line = br.readLine();
			if(line==null || !line.startsWith("#")) throw new IOException("Cannot read first line or "+uri);
			line=line.substring(1);//remove '#'
			final FileHeader header = new FileHeader(CharSplitter.TAB.split(line));
			while((line=br.readLine())!=null) {
				final Map<String,String> rec = header.toMap(CharSplitter.TAB.split(line));
				if(!rec.get("tax_id").equals("9606")) continue;
				final NcbiGeneInfo info =this.geneid2gene.getOrDefault(rec.get("GeneID"),null);
				if(info==null) {
					continue;
					}
				gene2go.add(new AbstractMap.SimpleEntry<>(info, rec.get("GO_ID")));
				}
			}
		final Comparator<Map.Entry<NcbiGeneInfo,String>> cmp = (A,B)->{
			return A.getKey().geneid.compareTo(B.getKey().geneid);
			};
		Collections.sort(gene2go,cmp);
		
		try(EqualRangeIterator<Map.Entry<NcbiGeneInfo,String>> eq = new EqualRangeIterator<>(gene2go.iterator(),cmp)) {
			while(eq.hasNext()) {
				final List<Map.Entry<NcbiGeneInfo,String>> array = eq.next();
				
				writer.writeStartElement("rdf", "Description", RDF.getURI());
				writer.writeAttribute("rdf", RDF.getURI(),"about",array.get(0).getKey().getURI());

				for(Map.Entry<NcbiGeneInfo,String> kv: array) {
					writer.writeEmptyElement(PREFIX, "has_go", NS);
					writer.writeAttribute("rdf",RDF.getURI(),"resource",
							"http://purl.obolibrary.org/obo/" + kv.getValue().replace(':', '_'));
					}
				writer.writeEndElement();
				writer.writeCharacters("\n");
				}
			}
		}

	private Map.Entry<String, NcbiGeneInfo> parseBiogridInteractor(StartElement root, final XMLEventReader xr) throws IOException,XMLStreamException {
		final QName QNAME_ID= new QName("id");
		Attribute att =root.getAttributeByName(QNAME_ID);
		String interactor_id = att.getValue();
		String ncbiTaxId = null;
		NcbiGeneInfo gene = null;
		while(xr.hasNext()) {
			final XMLEvent evt  = xr.nextEvent();
			if(evt.isStartElement()) {
				final StartElement startE = evt.asStartElement();
				final String lclName = startE.getName().getLocalPart();
				if(lclName.equals("organism")) {
					att = startE.getAttributeByName(new QName("ncbiTaxId"));
					if(att!=null) ncbiTaxId=att.getValue();
					}
				else if(lclName.equals("secondaryRef")) {
					att = startE.getAttributeByName(new QName("db"));
					if(att!=null && att.getValue().equals("entrez gene/locuslink")) {
						att = startE.getAttributeByName(QNAME_ID);
						if(att!=null) {
							final NcbiGeneInfo g = this.geneid2gene.get(att.getValue());
							if(g!=null) gene=g;
							}
						}
					}
				}
			else if(evt.isEndElement()) {
				final String lclName = evt.asEndElement().getName().getLocalPart();
				if(lclName.equals("proteinInteractor")) {
					break;
					}
				}
			}
		if(ncbiTaxId==null || !ncbiTaxId.equals("9606")) return null;
		if(gene==null) return null;
		return  new AbstractMap.SimpleEntry<String,NcbiGeneInfo>(interactor_id,gene);
		}
	
	private void parseBiogridInteraction(StartElement root, final XMLEventReader xr, final Map<String,NcbiGeneInfo> id2gene) throws IOException,XMLStreamException {
		Attribute att  = null;
		final Set<NcbiGeneInfo> genes = new HashSet<>();
		final Set<String> pmids = new HashSet<>();
		while(xr.hasNext()) {
			final XMLEvent evt  = xr.nextEvent();
			if(evt.isStartElement()) {
				final StartElement startE = evt.asStartElement();
				final String lclName = startE.getName().getLocalPart();
				if(lclName.equals("interactorRef")) {
					final String interactorRef = xr.getElementText();
					final NcbiGeneInfo gene = id2gene.get(interactorRef);
					if(gene!=null) genes.add(gene);
					}
				else if(lclName.equals("primaryRef")) {
					att = startE.getAttributeByName(new QName("db"));
					if(att!=null && att.getValue().equals("pubmed")) {
						att = startE.getAttributeByName(new QName("id"));
						if(att!=null) {
							pmids.add(att.getValue());
							}
						}
					}
				}
			else if(evt.isEndElement()) {
				final String lclName = evt.asEndElement().getName().getLocalPart();
				if(lclName.equals("interaction")) {
					break;
					}
				}
			}
		if(genes.size()>1) {
			writer.writeStartElement(PREFIX,"Interaction",NS);
			for(String pmid: pmids) {
				writer.writeEmptyElement(PREFIX, "pubmed", NS);
				writer.writeAttribute("rdf",RDF.getURI(),
						"resource",
						"http://www.ncbi.nlm.nih.gov/pubmed/"+pmid
						);
				}
			for(NcbiGeneInfo gene: genes) {
				writer.writeEmptyElement(PREFIX, "gene", NS);
				writer.writeAttribute("rdf",RDF.getURI(),
						"resource",
						gene.getURI()
						);
				
				}
			writer.writeEndElement();
			writer.writeCharacters("\n");
			}
		}
	
	private void parseBioGrid(String uri) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uri)) return;
		Path tmpZip = null;
		try {
			tmpZip = Files.createTempFile("biogrid", ".zip");
			try(InputStream in= IOUtils.openURIForReading(uri)) {
				IOUtils.copyTo(in, tmpZip);
				}
			final Map<String,NcbiGeneInfo> id2gene =new HashMap<>();
			try(InputStream in = Files.newInputStream(tmpZip)) {
				ZipInputStream zin = new ZipInputStream(in);
				for(;;) {
					final ZipEntry entry = zin.getNextEntry();
					if(entry==null) break;
					if(!entry.getName().endsWith(".psi.xml")) continue;
					final XMLInputFactory xif = XMLInputFactory.newFactory();
					final XMLEventReader xr = xif.createXMLEventReader(zin, "UTF-8");
					while(xr.hasNext()) {
						final XMLEvent evt  = xr.nextEvent();
						if(evt.isStartElement()) {
							final String lclName = evt.asStartElement().getName().getLocalPart();
							if(lclName.equals("proteinInteractor")) {
								final Map.Entry<String, NcbiGeneInfo> pair = parseBiogridInteractor(evt.asStartElement(), xr);
								if(pair!=null) {
									id2gene.put(pair.getKey(), pair.getValue());
									}
								}
							else if(lclName.equals("interaction")) {
								parseBiogridInteraction(evt.asStartElement(), xr, id2gene);							
								}

							}
						}
					xr.close();
					break;
					}
				}
			}
		catch(Throwable err) {
			if(tmpZip!=null) Files.delete(tmpZip);
			}
		}
	
	private void parseGFF(String uri,String build) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uri)) return;
		LOG.info("parsing "+uri);
		final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			final LineIterator li = LineIterators.of(br);
			codec.readHeader(li);
			while(!codec.isDone(li)) {
				Gff3Feature feat = codec.decode(li);
				if(feat==null || !feat.getType().equals("gene")) continue;
				final NcbiGeneInfo info = feat.getAttribute("gene_name").stream().
						map(F->this.symbol2gene.get(F)).
						filter(G->G!=null).
						findAny().
						orElse(null);
				if(info==null) {
					continue;
					}
				
				
				
				writer.writeStartElement("rdf", "Description", RDF.getURI());
				writer.writeAttribute("rdf", RDF.getURI(),"about",info.getURI());
				
				String geneid = feat.getAttribute("gene_id").stream().
					findAny().
					orElse(null);
				
				if(!StringUtils.isBlank(geneid) && geneid.startsWith("ENSG")) {
					final int dot = geneid.lastIndexOf(".");
					if(dot!=-1) geneid= geneid.substring(0,dot);
					writer.writeStartElement(PREFIX,"ensembl_geneid",NS);
					writer.writeCharacters(geneid);
					writer.writeEndElement();
					}
				
				String hgnc_id = feat.getAttribute("hgnc_id").stream().
						findAny().
						orElse(null);
				if(!StringUtils.isBlank(hgnc_id)) {
					writer.writeStartElement(PREFIX,"hgnc_id",NS);
					writer.writeCharacters(hgnc_id);
					writer.writeEndElement();
				}
				
				writer.writeStartElement(PREFIX,"location",NS);
				
				writer.writeStartElement(PREFIX,"Interval",NS);
				
				writer.writeStartElement(PREFIX,"build",NS);
				writer.writeCharacters(build);
				writer.writeEndElement();

				
				writer.writeStartElement(PREFIX,"contig",NS);
				writer.writeCharacters(feat.getContig());
				writer.writeEndElement();
				
				writer.writeStartElement(PREFIX,"start",NS);
				writer.writeAttribute("rdf", RDF.getURI(), "datatype", XSD.integer.getURI());
				writer.writeCharacters(String.valueOf(feat.getStart()-1));
				writer.writeEndElement();

				writer.writeStartElement(PREFIX,"end",NS);
				writer.writeAttribute("rdf", RDF.getURI(), "datatype",XSD.integer.getURI());
				writer.writeCharacters(String.valueOf(feat.getEnd()));
				writer.writeEndElement();

				
				writer.writeEndElement(); // Interval
				writer.writeEndElement();//location
				writer.writeEndElement();//resource
				writer.writeCharacters("\n");
				}
			codec.close(li);	
		}
	}
	
	private EventFilter createOWLEventFilter() {
		return  new SimpleEventFilter(STACK->{
			final String lcl = STACK.get(STACK.size()-1).getLocalPart();
			if(lcl.equals("Ontology")) return false;
			if(lcl.equals("AnnotationProperty")) return false;
			if(lcl.equals("Axiom")) return false;
			//if(lcl.equals("hasDbXref")) return false;
			if(lcl.equals("hasRelatedSynonym")) return false;
			if(lcl.equals("hasExactSynonym")) return false;
			if(lcl.equals("hasOBONamespace")) return false;
			if(lcl.equals("inSubset")) return false;
		return true;
		});
	}
	
	
	private void parseOBOOWLOntology(
			final String uri,
			final String className,
			final Predicate<String> acceptTermId
			) throws IOException,XMLStreamException,SAXException {

		if(StringUtils.isBlank(uri)) {
			LOG.info("skiping OWL ontology... for "+PREFIX+":"+className+" because URL is blank.");
			return;
			}
		
		this.writer.writeComment("BEGIN parsing ontology "+uri);
		this.writer.writeCharacters("\n");
		
		final Model ontModel  = ModelFactory.createDefaultModel();
		final XMLInputFactory inputFactory = XMLInputFactory.newInstance();
		inputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
		
		try(Reader r= IOUtils.openURIForBufferedReading(uri)) {
			final XMLEventReader reader0 = inputFactory.createXMLEventReader(r);
			final XMLEventReader reader = inputFactory.createFilteredReader(reader0,createOWLEventFilter());
		    StAX2Model.read(reader,ontModel,IOUtil.isUrl(uri)?uri:"file://"+uri); 
			reader0.close();
			}
		if(ontModel.size()==0) {
			LOG.warn("no statement foudn in "+uri);
			return ;
			}
		final ExtendedIterator<Resource> r1= ontModel.listSubjectsWithProperty(RDF.type, OWL.Class);
		final Property deprecated = ResourceFactory.createProperty(OWL.getURI(), "deprecated");
		while(r1.hasNext()) {
			final Resource rsrc = r1.next();
			if(rsrc.hasLiteral(deprecated, true)) continue;
			final String term_id = JenaUtils.stream(rsrc.listProperties(OBOInOwl.id)).
					map(R->R.getObject()).
					filter(R->R.isLiteral()).
					map(R->R.asLiteral()).
					map(L->L.getString()).
					findFirst().
					orElse(null);
			if(term_id==null) continue;
			if(!acceptTermId.test(term_id)) continue;
			
			
			writer.writeStartElement(PREFIX,className ,NS);
			writer.writeAttribute("rdf", RDF.getURI(), "about", rsrc.getURI() );
			
			String s = JenaUtils.stream(rsrc.listProperties(RDFS.label)).
					map(R->R.getObject()).
					filter(R->R.isLiteral()).
					map(R->R.asLiteral()).
					map(L->L.getString()).
					findFirst().
					orElse(null);
			
			if(!StringUtils.isBlank(s)) {
				writer.writeStartElement("dc","description",DC.getURI());
				writer.writeCharacters(s);
				writer.writeEndElement();
				}
		
			
			
			if(!StringUtils.isBlank(s)) {
				writer.writeStartElement("dc","title",DC.getURI());
				writer.writeCharacters(term_id);
				writer.writeEndElement();
				}
			
			
			final ExtendedIterator<Statement> r2= rsrc.listProperties(RDFS.subClassOf);
			while(r2.hasNext()) {
				final Statement stmt = r2.next();
				if(!stmt.getObject().isResource()) continue;
				if(stmt.getObject().isAnon()) continue;
				if(stmt.getObject().equals(rsrc)) continue;
				writer.writeEmptyElement("rdfs", "subClassOf", RDFS.getURI());
				writer.writeAttribute("rdf",RDF.getURI(),
						"resource",
						stmt.getObject().asResource().toString());
				}
			r2.close();
			
			writer.writeEndElement();
			writer.writeCharacters("\n");
			}
		r1.close();
		
		this.writer.writeComment("END parsing ontology "+uri);
		this.writer.writeCharacters("\n");

		}
	
	private void parseHumanDiseaseOntology(String uri) throws IOException,XMLStreamException,SAXException {
		parseOBOOWLOntology(uri, "Disease" ,S->true);
		}
	
	private void parseHumanPhenotypeOntology(String uri) throws IOException,XMLStreamException,SAXException {
		parseOBOOWLOntology(uri, "Phenotype",S->S.startsWith("HP:"));
		}
	private void parseGO(String uri) throws IOException,XMLStreamException,SAXException {
		parseOBOOWLOntology(uri, "GOTerm",S->S.startsWith("GO:"));
		}

	private void parseMondoOWL(String uri)  throws IOException,XMLStreamException,SAXException {
		parseOBOOWLOntology(uri, "Disease" ,S->true);
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			if(!args.isEmpty()) {
				LOG.error("Too many arguments.");
				return -1;
				}
			
		
			final String encoding="UTF-8";
			final XMLOutputFactory xof = XMLOutputFactory.newFactory();
			try(BufferedWriter bw = (outputFile==null?
					new BufferedWriter(new OutputStreamWriter(stdout(), encoding)) :
					Files.newBufferedWriter(this.outputFile,Charset.forName(encoding))
					)) {
				this.writer = xof.createXMLStreamWriter(bw);
				this.writer.writeStartDocument(encoding, "1.0");
				this.writer.writeStartElement("rdf", "RDF", RDF.getURI());
				this.writer.writeNamespace(PREFIX, NS);
				this.writer.writeNamespace("rdf",RDF.getURI());
				this.writer.writeNamespace("rdfs",RDFS.getURI());
				this.writer.writeNamespace("xsd",XSD.getURI());
				this.writer.writeNamespace("dc",DC.getURI());
				this.writer.writeCharacters("\n");
				this.writer.flush();
				
				//parseMondoOWL(resourceMap.getOrDefault("MONDO_OWL",""));
				
				/** NCBI gene INFO */
				this.writer.writeComment("NCBI GENE INFO");
				this.writer.writeCharacters("\n");
				
				parseNcbiGeneInfo(resourceMap.getOrDefault("NCBI_GENE_INFO",""));
				
				parseHumanPhenotypeOntology(resourceMap.getOrDefault("HUMAN_HPO_OWL", ""));
				parseHPO2gene(resourceMap.getOrDefault("HPO_PHENOTYPE_TO_GENE", ""));
				parseNcbiGeneGO(resourceMap.getOrDefault("NCBI_GENE_GO",""));
				parseGO(resourceMap.getOrDefault("GO_OWL",""));
				
				final String gencode_release = resourceMap.getOrDefault("GENCODE_RELEASE", "43");
				if(!StringUtils.isBlank(gencode_release)) {
					parseGFF(resourceMap.getOrDefault("GFF3_GRCH38","").replace("{GENCODE_RELEASE}", gencode_release),
							"grch38"
							);
					parseGFF(resourceMap.getOrDefault("GFF3_GRCH37","").replace("{GENCODE_RELEASE}", gencode_release),
							"grch37"
							);
					}
				
				parseBioGrid(resourceMap.getOrDefault("BIOGRID_XML_25",""));

				//parseHumanDiseaseOntology(resourceMap.getOrDefault("HUMAN_DO_OWL",""));
				
				
				this.writer.writeEndElement();//RDF
				this.writer.flush();
				this.writer.close();
				this.writer = null;
				bw.flush();
				}
			return 0;
			}
		catch (Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(final String[] args) {
		new BioToRDF().instanceMainWithExit(args);
	}
	
}
