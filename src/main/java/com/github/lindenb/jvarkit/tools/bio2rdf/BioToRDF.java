package com.github.lindenb.jvarkit.tools.bio2rdf;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

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
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.iterator.LineIterators;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.xml.SimpleEventFilter;

import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.readers.LineIterator;

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
	@DynamicParameter(names={"-D"},description="parameters. -Dkey1=value1  -Dkey2=value2 ...")
	private Map<String,String> dynaParams = new HashMap<String,String>() {{{
		put("NCBI_GENE_INFO", "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz");
		put("NCBI_GENE_GO","https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz");
		put("GENCODE_RELEASE", "43");
		put("GFF3_GRCH38", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.annotation.gff3.gz");
		put("GFF3_GRCH37", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/GRCh37_mapping/gencode.v{GENCODE_RELEASE}lift37.annotation.gff3.gz");
		put("HUMAN_DO_OWL", "https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/main/src/ontology/HumanDO.owl");
		}}};


	private XMLStreamWriter writer = null;
	private final Map<String, NcbiGeneInfo> geneid2gene = new HashMap<>();
	private final Map<String, NcbiGeneInfo> symbol2gene = new HashMap<>();
	
	private static class NcbiGeneInfo {
		String geneid;
		String symbol;
		public String getURI() {
			return "https://www.ncbi.nlm.nih.gov/gene/"+ this.geneid;
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
		try(BufferedReader br = IOUtils.openURIForBufferedReading(uri)) {
			String line = br.readLine();
			if(line==null || !line.startsWith("#")) throw new IOException("Cannot read first line or "+uri);
			line=line.substring(1);//remove '#'
			final FileHeader header = new FileHeader(CharSplitter.TAB.split(line));
			while((line=br.readLine())!=null) {
				final Map<String,String> rec = header.toMap(CharSplitter.TAB.split(line));
				if(!rec.get("tax_id").equals("9606")) continue;
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
				if(!StringUtils.isBlank(s)) {
					writer.writeStartElement(PREFIX,"gene_type",NS);
					writer.writeCharacters(s.replace('-', '_'));
					writer.writeEndElement();
					}
				
 				for(String syn: CharSplitter.PIPE.split(rec.get("Synonyms"))) {
 					if(StringUtils.isBlank(syn)) continue;
 					writer.writeStartElement(PREFIX,"synonym",NS);
					writer.writeCharacters(syn);
					writer.writeEndElement();
 					}
 				for(String xref: CharSplitter.PIPE.split(rec.get("dbXrefs"))) {
 					if(StringUtils.isBlank(xref)) continue;
 					int colon = xref.indexOf(':');
 					if(colon==-1) continue;
 					String key = xref.substring(0,colon);
 					if(key.equals("Ensembl")) {
 						writer.writeStartElement(PREFIX,"ensembl_geneid",NS);
 						writer.writeCharacters(xref.substring(colon+1));
 						writer.writeEndElement();
 						}
 					else if(key.equals("HGNC")) {
 						writer.writeStartElement(PREFIX,"hgnc_id",NS);
 						writer.writeCharacters(xref.substring(colon+1));
 						writer.writeEndElement();
 						}
 					}
				writer.writeEndElement();
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
				
				writer.writeStartElement("rdf", "Resource", RDF.getURI());
				writer.writeAttribute("rdf", RDF.getURI(),"about",array.get(0).getKey().getURI());

				for(Map.Entry<NcbiGeneInfo,String> kv: array) {
					writer.writeEmptyElement(PREFIX, "has_go", NS);
					writer.writeAttribute("rdf",RDF.getURI(),"resource",
							"http://purl.obolibrary.org/obo/" + kv.getValue().replace(':', '_'));
					}
				writer.writeEndElement();
				}
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
				
				
				
				writer.writeStartElement("rdf", "Resource", RDF.getURI());
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
				
				writer.writeStartElement(PREFIX,"Location",NS);
				
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

				
				writer.writeEndElement(); // Location
				writer.writeEndElement();//location
				writer.writeEndElement();//resource
				}
			codec.close(li);	
		}
	}
	
	private void parseHumanDiseaseOntology(String uri) throws IOException,XMLStreamException,SAXException {
		if(StringUtils.isBlank(uri)) return;
		final Model ontModel  = ModelFactory.createDefaultModel();
		final XMLInputFactory inputFactory = XMLInputFactory.newInstance();
		inputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
		try(Reader r= IOUtils.openURIForBufferedReading(uri)) {
			final XMLEventReader reader0 = inputFactory.createXMLEventReader(r);
			final XMLEventReader reader = inputFactory.createFilteredReader(reader0, new SimpleEventFilter(qName->{
					final String lcl = qName.getLocalPart();
					if(lcl.equals("Ontology")) return false;
					if(lcl.equals("AnnotationProperty")) return false;
					if(lcl.equals("Axiom")) return false;
					if(lcl.equals("hasDbXref")) return false;
					if(lcl.equals("hasExactSynonym")) return false;
					if(lcl.equals("hasOBONamespace")) return false;
					if(lcl.equals("inSubset")) return false;
				return true;
				}));
		    StAX2Model.read(reader,ontModel,"file://"+uri); 
			reader0.close();
			}
		final String OBO_IN_OWL= "http://www.geneontology.org/formats/oboInOwl#";
		ExtendedIterator<Resource> r1= ontModel.listSubjectsWithProperty(RDF.type, OWL.Class);
		final Property deprecated = ResourceFactory.createProperty(OWL.getURI(), "deprecated");
		final Property oboid = ResourceFactory.createProperty(OBO_IN_OWL, "id");
		while(r1.hasNext()) {
			final Resource rsrc = r1.next();
			if(rsrc.hasLiteral(deprecated, true)) continue;
			writer.writeStartElement(PREFIX,"Disease",NS);
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
		
			s = JenaUtils.stream(rsrc.listProperties(oboid)).
					map(R->R.getObject()).
					filter(R->R.isLiteral()).
					map(R->R.asLiteral()).
					map(L->L.getString()).
					findFirst().
					orElse(null);
			
			if(!StringUtils.isBlank(s)) {
				writer.writeStartElement("dc","title",DC.getURI());
				writer.writeCharacters(s);
				writer.writeEndElement();
				}
			
			
			final ExtendedIterator<Statement> r2= rsrc.listProperties(RDFS.subClassOf);
			while(r2.hasNext()) {
				final Statement rsrc2 = r2.next();
				if(!rsrc2.getObject().isResource()) continue;
				
				writer.writeEmptyElement("rdfs", "subClassOf", RDFS.getURI());
				writer.writeAttribute("rdf",RDF.getURI(),"resource",rsrc2.getObject().asResource().getURI());
				
				}
			r2.close();
			
			writer.writeEndElement();
			writer.writeCharacters("\n");
			}
		r1.close();
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final Map<String,String> resourceMap = new HashMap<>();
			try(BufferedReader br = super.openBufferedReader(oneFileOrNull(args))) {
				String line;
				while((line=br.readLine())!=null) {
					if(StringUtils.isBlank(line)) continue;
					final int tab = line.indexOf('\t');
					if(tab==-1) {
						LOG.error("tab missing in resources file: "+line);
						return -1;
						}
					final String key = line.substring(0,tab).toUpperCase().trim();
					if(resourceMap.containsKey(key)) {
						LOG.error("duplicate resource key "+line);
						return -1;
						}
					resourceMap.put(
							key,
							line.substring(tab).trim()
							);
					}
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
				this.writer.writeCharacters("\n");
				
				/** NCBI gene INFO */
				this.writer.writeComment("NCBI GENE INFO");
				this.writer.writeCharacters("\n");
				
				parseNcbiGeneInfo(resourceMap.getOrDefault("NCBI_GENE_INFO",""));
				parseNcbiGeneGO(resourceMap.getOrDefault("NCBI_GENE_GO",""));
				
				final String gencode_release = resourceMap.getOrDefault("GENCODE_RELEASE", "43");
				if(!StringUtils.isBlank(gencode_release)) {
					parseGFF(resourceMap.getOrDefault("GFF3_GRCH38","").replace("{GENCODE_RELEASE}", gencode_release),
							"grch38"
							);
					parseGFF(resourceMap.getOrDefault("GFF3_GRCH37","").replace("{GENCODE_RELEASE}", gencode_release),
							"grch37"
							);
					}
				
				
				//parseHPOA(resourceMap.getAttribute("HPOA", "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"));
				parseHumanDiseaseOntology(resourceMap.getOrDefault("HUMAN_DO_OWL",""));
				
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
	
	public static void main(String[] args) {
		new BioToRDF().instanceMainWithExit(args);
	}
	
}
