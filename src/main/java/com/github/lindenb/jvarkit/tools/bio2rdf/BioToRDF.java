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
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import javax.xml.stream.EventFilter;
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
import com.github.lindenb.jvarkit.jena.vocabulary.OBOInOwl;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
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


END_DOC
 
**/
@Program(
		name="bio2rdf",
		description="Build a RDF database for human from misc sources",
		keywords={"rdf","ontology","sparql"},
		creationDate="20220427",
		modificationDate="20220510",
		jvarkit_amalgamion =  true
		)
public class BioToRDF extends Launcher {
	private static final Logger LOG = Logger.of(BioToRDF.class);
	private final static String PREFIX="bio";
	private final static String NS = "https://umr1087.univ-nantes.fr/bio2rdf/";
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"--genes"},description="Limit to those genes names , separated with comma (for debugging)")
	private String limitGenesStr="";
	@Parameter(names={"--min-stringdb-combined-score"},description="Discard interaction with stringdb combined score < 'x' ")
	private int min_stringdb_combined_score=990;

	
	
	
	@SuppressWarnings("serial")
	@DynamicParameter(names={"-D"},description="parameters. -Dkey1=value1  -Dkey2=value2 ...")
	private Map<String,String> resourceMap = new HashMap<String,String>() {{{
		put("NCBI_GENE_INFO", "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz");
		put("NCBI_GENE_GO","https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz");
		put("GENCODE_RELEASE", "43");
		put("GFF3_GRCH38", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/gencode.v{GENCODE_RELEASE}.annotation.gff3.gz");
		put("GFF3_GRCH37", "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{GENCODE_RELEASE}/GRCh37_mapping/gencode.v{GENCODE_RELEASE}lift37.annotation.gff3.gz");
		//put("HUMAN_DO_OWL", "https://github.com/DiseaseOntology/HumanDiseaseOntology/raw/main/src/ontology/HumanDO.owl");
		//put("HUMAN_HPO_OWL", "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/hp.owl");
		put("HPO_PHENOTYPE_TO_GENE", "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2023-04-05/phenotype_to_genes.txt");
		put("GO_OWL","http://purl.obolibrary.org/obo/go.owl");
		put("STRINGDB_RELEASE","11.5");
		put("STRINGDB_PROTEIN_ALIASES","https://stringdb-static.org/download/protein.aliases.v{STRINGDB_RELEASE}/9606.protein.aliases.v{STRINGDB_RELEASE}.txt.gz");
		put("STRINGDB_LINK","https://stringdb-static.org/download/protein.links.v{STRINGDB_RELEASE}/9606.protein.links.v{STRINGDB_RELEASE}.txt.gz");
		//put("MONDO_OWL","https://github.com/monarch-initiative/mondo/releases/download/v2023-04-04/mondo.owl");
		put("BASE_NCBI_REFSEQ","https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene");
		}}};


	private XMLStreamWriter writer = null;
	private final Map<String, NcbiGeneInfo> geneid2gene = new HashMap<>();
	private final Map<String, NcbiGeneInfo> symbol2gene = new HashMap<>();
	private final Map<String, NcbiGeneInfo> hgnc2gene = new HashMap<>();
	private static class Build {
		final String uri;
		final String label;
		Build(final String uri,final String label) {
			this.uri = uri;
			this.label = label;
			}
		}
	
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
	
	private String log(String uri) {
		LOG.info("Parsing "+uri);
		return uri;
		}
	
	
	private void parseNcbiGeneInfo(String uri) throws IOException,XMLStreamException {
		
		final Set<String> limit_gene_names = Arrays.stream(limitGenesStr.split("[, \t]")).
				filter(S->!StringUtils.isBlank(S)).
				collect(Collectors.toSet());
				
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(uri))) {
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

				for(String xref: CharSplitter.PIPE.split(rec.get("dbXrefs"))) {
 					if(StringUtils.isBlank(xref)) continue;
 					final int colon = xref.indexOf(':');
 					if(colon==-1) continue;
 					final String key = xref.substring(0,colon);
 					if(key.equals("HGNC")) {
 						final String hgnc_id= xref.substring(colon+1);
 						info.hgnc = hgnc_id;
 						break;
 						}
 					}
				
				this.geneid2gene.put(info.geneid, info);
				this.symbol2gene.put(info.symbol, info);
				
				writer.writeStartElement("rdf","Resource",RDF.getURI());
				writer.writeAttribute("rdf", RDF.getURI(), "about", info.getURI() );
				
				writer.writeStartElement("dc","title",DC.getURI());
				writer.writeCharacters(info.symbol);
				writer.writeEndElement();
				
				writer.writeStartElement(PREFIX,"symbol",NS);
				writer.writeCharacters(info.symbol);
				writer.writeEndElement();
				writer.writeStartElement(PREFIX,"ncbi_gene_id",NS);
				writer.writeCharacters(info.geneid);
				writer.writeEndElement();

				if(!StringUtils.isBlank(info.hgnc)) {
					writer.writeStartElement(PREFIX,"hgnc_id",NS);
					writer.writeCharacters(info.hgnc);
					writer.writeEndElement();
					hgnc2gene.put(info.hgnc, info);
					}
				
				
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
 					}
				writer.writeEndElement();
				writer.writeCharacters("\n");
				}
			}
		}
	
	
	private void parseNcbiRefseq(final String baserefseq) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(baserefseq)) return;
		final List<String> gb_uris = new ArrayList<>();
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(baserefseq + "/refseqgene.files.installed"))) {
			br.lines().
				map(S-> CharSplitter.TAB.split(S)).
				filter(T->T[1].endsWith(".gbff.gz")).
				forEach(T->gb_uris.add(baserefseq+"/"+T[1]));
			}
		
		if(gb_uris.isEmpty()) throw new IOException("cannot any url in "+baserefseq + "/refseqgene.files.installed");
		
		final Map<String,NcbiGeneInfo> refseq2gene= new HashMap<>();
		final String LRG_RefSeqGene = baserefseq+"/LRG_RefSeqGene";
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(LRG_RefSeqGene))) {
			String line=br.readLine();
			if(line==null) throw new IOException("cannot get first line of "+LRG_RefSeqGene);
			final FileHeader header = new FileHeader(CharSplitter.TAB.split(line));
			while((line=br.readLine())!=null) {
				final Map<String,String> rec = header.toMap(CharSplitter.TAB.split(line));
				String acn = rec.get("RSG");
				int dot = acn.lastIndexOf('.');
				if(dot!=-1) acn=acn.substring(0,dot);
				final NcbiGeneInfo info = this.geneid2gene.get(rec.get("GeneID"));

				if(info!=null) refseq2gene.put(acn, info);
				}
			}
		
		for(final String uri :gb_uris )  {
			if(refseq2gene.isEmpty()) break;
			try(BufferedReader br = IOUtils.openURIForBufferedReading(log(uri))) {
				String line;
				String locus="";
				int state=-1;
				StringBuilder desc=new StringBuilder();
				while((line=br.readLine())!=null) {
					if(line.startsWith("//")) {
						final NcbiGeneInfo info = refseq2gene.get(locus);
						if(info!=null) {
							writer.writeStartElement("rdf", "Description", RDF.getURI());
							writer.writeAttribute("rdf", RDF.getURI(),"about",info.getURI());

							writer.writeStartElement(PREFIX, "refseq_summary", NS);
							writer.writeCharacters(StringUtils.normalizeSpaces(desc.toString()));
							writer.writeEndElement();
							
							writer.writeEndElement();
							writer.writeCharacters("\n");
							
							refseq2gene.remove(locus);//faster ?
							}
						state=-1;
						locus="";
						desc.setLength(0);
						}
					else if(line.startsWith("VERSION ")) {
						locus = line.substring(12);
						int dot = locus.lastIndexOf('.');
						if(dot!=-1) locus=locus.substring(0,dot);
						state=-1;
						}
					else if(line.startsWith("COMMENT ")) {
						desc.append(line.substring(12));
						state = 1;
						}
					else if(line.startsWith("  ") && state==1) {
						desc.append(line.substring(12));
						}
					else {
						state=-1;
						}
					}
				}
			}
		for(final String key: refseq2gene.keySet()) {
			LOG.warn("cannot find refseq for "+key);
			}
		}
	
	private void parseNcbiGeneGO(final String uri) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uri)) return;
		final List<Map.Entry<NcbiGeneInfo,String>> gene2go = new ArrayList<>();
		
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(uri))) {
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


	private void parseStringDB(String uriAliases, final String uriLinks) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uriAliases)) return;
		if(StringUtils.isBlank(uriLinks)) return;
		final Map<String,NcbiGeneInfo> proteinId2gene = new HashMap<>(20_000);
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(uriAliases))) {
			String line;
			while((line=br.readLine())!=null) {
				final String[] tokens = CharSplitter.TAB.split(line);
				if(tokens.length!=3) throw new JvarkitException.TokenErrors(3, tokens);
				if(!tokens[2].equals("Ensembl_HGNC_HGNC_ID")) continue;
				final NcbiGeneInfo gene = this.hgnc2gene.get(tokens[1]);
				if(gene==null) continue;
				proteinId2gene.put(tokens[0], gene);
				}
			}
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(uriLinks))) {
			String line = br.readLine();
			if(line==null) throw new IOException("Cannot read first line of "+uriLinks);
			final FileHeader header = new FileHeader(CharSplitter.SPACE.split(line));
			while((line=br.readLine())!=null) {
				final Map<String,String> row = header.toMap(CharSplitter.SPACE.split(line));
				
				final int combined_score = Integer.parseInt(row.get("combined_score"));
				if(combined_score < this.min_stringdb_combined_score) continue;
				
				
				final NcbiGeneInfo gene1 = proteinId2gene.get(row.get("protein1"));
				if(gene1==null) continue;
				final NcbiGeneInfo gene2 = proteinId2gene.get(row.get("protein2"));
				if(gene2==null) continue;
				writer.writeStartElement(PREFIX,"StringDBInteraction",NS);
				
				for(NcbiGeneInfo gene:Arrays.asList(gene1,gene2)) {
					writer.writeEmptyElement(PREFIX, "interactor", NS);
					writer.writeAttribute("rdf",RDF.getURI(),"resource", gene.getURI());
					}
				
				writer.writeStartElement(PREFIX,"score",NS);
				writer.writeAttribute("rdf", RDF.getURI(), "datatype", XSD.integer.getURI());
				writer.writeCharacters(row.get("combined_score"));
				writer.writeEndElement();

				
				writer.writeComment(gene1.symbol+" " + gene2.symbol);
				
				writer.writeEndElement();
				writer.writeCharacters("\n");
				}
			}
		}
	
	private void parseGFF(final String uri,final Build build) throws IOException,XMLStreamException {
		if(StringUtils.isBlank(uri)) return;
		LOG.info("parsing "+uri);
		
		writer.writeComment("PROCESSING GFF "+uri);
		
		writer.writeStartElement(PREFIX, "Build", NS);
		writer.writeAttribute("rdf", RDF.getURI(),"about",build.uri);
		writer.writeStartElement("rdfs","label",RDFS.getURI());
		writer.writeCharacters(build.label);
		writer.writeEndElement();
		writer.writeEndElement();
		writer.writeCharacters("\n");
		
		final Gff3Codec codec = new Gff3Codec(Gff3Codec.DecodeDepth.SHALLOW);
		try(BufferedReader br = IOUtils.openURIForBufferedReading(log(uri))) {
			final LineIterator li = LineIterators.of(br);
			codec.readHeader(li);
			while(!codec.isDone(li)) {
				final Gff3Feature feat = codec.decode(li);
				if(feat==null || !feat.getType().equals("gene")) continue;
				NcbiGeneInfo info =  feat.getAttribute("hgnc_id").stream().
						map(F->this.hgnc2gene.get(F)).
						filter(G->G!=null).
						findAny().
						orElse(null);
				if(info==null) {
						info = feat.getAttribute("gene_name").stream().
							map(F->this.symbol2gene.get(F)).
							filter(G->G!=null).
							findAny().
							orElse(null);
					}
				
				
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
				
				writer.writeEmptyElement(PREFIX,"build",NS);
				writer.writeAttribute("rdf", RDF.getURI(), "resource",build.uri);


				
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
		
		try(Reader r= IOUtils.openURIForBufferedReading(log(uri))) {
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
	
	
	
	private void parseGO(String uri) throws IOException,XMLStreamException,SAXException {
		parseOBOOWLOntology(uri, "GOTerm",S->S.startsWith("GO:"));
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
				
				parseNcbiRefseq(resourceMap.getOrDefault("BASE_NCBI_REFSEQ", ""));
				
				final String string_db_release = resourceMap.getOrDefault("STRINGDB_RELEASE", "11.5");;
				parseStringDB(
						resourceMap.getOrDefault("STRINGDB_PROTEIN_ALIASES","").replace("{STRINGDB_RELEASE}", string_db_release),
						resourceMap.getOrDefault("STRINGDB_LINK","").replace("{STRINGDB_RELEASE}", string_db_release)
						);

				
				parseNcbiGeneGO(resourceMap.getOrDefault("NCBI_GENE_GO",""));
				parseGO(resourceMap.getOrDefault("GO_OWL",""));
				
				final String gencode_release = resourceMap.getOrDefault("GENCODE_RELEASE", "43");
				if(!StringUtils.isBlank(gencode_release)) {
					parseGFF(resourceMap.getOrDefault("GFF3_GRCH38","").replace("{GENCODE_RELEASE}", gencode_release),
							new Build("https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.40","grch38")
							);
					parseGFF(resourceMap.getOrDefault("GFF3_GRCH37","").replace("{GENCODE_RELEASE}", gencode_release),
							new Build("https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13","grch37")
							);
					}
				
				
				
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
