/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.go;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gexf.GexfConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.go.GOOntology;
import com.github.lindenb.jvarkit.go.GOParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;
import com.github.lindenb.jvarkit.goa.GOAFileIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.gff.Gff3Codec;
import htsjdk.tribble.gff.Gff3Codec.DecodeDepth;
import htsjdk.tribble.gff.Gff3Feature;
import htsjdk.tribble.gff.Gff3Writer;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;

/**
BEGIN_DOC

## Example

children of  GO:0005216 'ion channel activity' 

```
$ java -jar dist/goutils.jar -R is_a  -A 'GO:0005216' 

#ACN	NAME	DEFINITION
GO:1905030	voltage-gated ion channel activity involved in regulation of postsynaptic membrane potential	Any voltage-gated ion channel activity that is involved in regulation of postsynaptic membrane potential.
GO:1905057	voltage-gated calcium channel activity involved in regulation of postsynaptic cytosolic calcium levels	Any voltage-gated calcium channel activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:1905054	calcium-induced calcium release activity involved in regulation of presynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of presynaptic cytosolic calcium ion concentration.
GO:1905058	calcium-induced calcium release activity involved in regulation of postsynaptic cytosolic calcium ion concentration	Any calcium-induced calcium release activity that is involved in regulation of postsynaptic cytosolic calcium ion concentration.
GO:0016286	small conductance calcium-activated potassium channel activity	Enables the transmembrane transfer of potassium by a channel with a unit conductance of 2 to 20 picoSiemens that opens in response to stimulus by internal calcium ions. Small conductance calcium-activated potassium channels are more sensitive to calcium than are large conductance calcium-activated potassium channels. Transport by a channel involves catalysis of facilitated diffusion of a solute (by an energy-independent process) involving passage through a transmembrane aqueous pore or channel, without evidence for a carrier-mediated mechanism.
GO:0043855	cyclic nucleotide-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0043854	cyclic nucleotide-gated mechanosensitive ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens in response to a mechanical stress and when a cyclic nucleotide has been bound by the channel complex or one of its constituent parts.
GO:0099142	intracellularly ATP-gated ion channel activity	Enables the transmembrane transfer of an ion by a channel that opens when ATP has been bound by the channel complex or one of its constituent parts on the intracellular side of the plasma membrane.
GO:0099101	G-protein gated potassium channel activity	A potassium channel activity that is gated by binding of a G-protein beta-gamma dimer.
(...)
```

## Example

### action = goa

```
$ wget -q -O - "http://geneontology.org/gene-associations/goa_human.gaf.gz" | gunzip -c | java -jar dist/goutils.jar -go go.obo --action goa -A 'ion transmembrane transport'  | head
UniProtKB	A0A1W2PN81	CHRNA7	enables	GO:0022848	GO_REF:0000002	IEA	InterPro:IPR002394	FNeuronal acetylcholine receptor subunit alpha-7	CHRNA7	protein	taxon:9606	20210612	InterPro
UniProtKB	A0PJK1	SLC5A10	enables	GO:0015370	Reactome:R-HSA-8876283	TAS		F	Sodium/glucose cotransporter 5	SLC5A10|SGLT5	protein	taxon:9606	20200515	Reactome
UniProtKB	A0PJK1	SLC5A10	involved_in	GO:0035725	GO_REF:0000108	IEA	GO:0015370	P	Sodium/glucose cotransporter 5	SLC5A10|SGLT5	protein	taxon:9606	20210613	GOC
UniProtKB	A1A4F0	SLC66A1L	involved_in	GO:1903401	GO_REF:0000108	IEA	GO:0015189	PPutative uncharacterized protein SLC66A1L	SLC66A1L|C3orf55|PQLC2L	protein	taxon:9606	20210613	GOC
UniProtKB	A1A4F0	SLC66A1L	involved_in	GO:1903826	GO_REF:0000108	IEA	GO:0015181	PPutative uncharacterized protein SLC66A1L	SLC66A1L|C3orf55|PQLC2L	protein	taxon:9606	20210613	GOC
UniProtKB	A1A5B4	ANO9	enables	GO:0005229	PMID:22946059	IMP		F	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20140424	UniProt
UniProtKB	A1A5B4	ANO9	enables	GO:0005229	Reactome:R-HSA-2684901	TAS		F	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20200515	Reactome
UniProtKB	A1A5B4	ANO9	involved_in	GO:0034220	Reactome:R-HSA-983712	TAS		P	Anoctamin-9	ANO9|PIG5|TMEM16J|TP53I5	protein	taxon:9606	20210310	Reactome
UniProtKB	A5X5Y0	HTR3E	enables	GO:0022850	PMID:17392525	IDA		F	5-hydroxytryptamine receptor 3E	HTR3E	protein	taxon:9606	20130210	CACAO
UniProtKB	A6NJY1	SLC9B1P1	enables	GO:0015299	GO_REF:0000002	IEA	InterPro:IPR006153	FPutative SLC9B1-like protein SLC9B1P1	SLC9B1P1	protein	taxon:9606	20210612	InterPro
```

### action = gff3

```
$ wget -q -O - "http://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/Homo_sapiens.GRCh37.82.chr.gff3.gz" | gunzip -c | java -jar dist/goutils.jar --action gff3 -A 'GO:0005216'   | head
WARNING	2021-10-21 16:19:29	AsciiLineReader	Creating an indexable source for an AsciiFeatureCodec using a stream that is neither a PositionalBufferedStream nor a BlockCompressedInputStream
##gff-version 3.1.25
1	ensembl_havana	gene	1215816	1227409	.	+	.	ID=gene%3AENSG00000162572;Name=SCNN1D;biotype=protein_coding;description=sodium channel%2C non-voltage-gated 1%2C delta subunit %5BSource%3AHGNC Symbol%3BAcc%3A10601%5D;gene_id=ENSG00000162572;logic_name=ensembl_havana_gene;version=15
1	ensembl_havana	gene	1950780	1962192	.	+	.	ID=gene%3AENSG00000187730;Name=GABRD;biotype=protein_coding;description=gamma-aminobutyric acid %28GABA%29 A receptor%2C delta %5BSource%3AHGNC Symbol%3BAcc%3A4084%5D;gene_id=ENSG00000187730;logic_name=ensembl_havana_gene;version=6
1	ensembl_havana	gene	6051526	6161253	.	+	.	ID=gene%3AENSG00000069424;Name=KCNAB2;biotype=protein_coding;description=potassium voltage-gated channel%2C shaker-related subfamily%2C beta member 2 %5BSource%3AHGNC Symbol%3BAcc%3A6229%5D;gene_id=ENSG00000069424;logic_name=ensembl_havana_gene;version=10
1	ensembl_havana	gene	11866207	11903201	.	+	.ID=gene%3AENSG00000011021;Name=CLCN6;biotype=protein_coding;description=chloride channel%2C voltage-sensitive 6 %5BSource%3AHGNC Symbol%3BAcc%3A2024%5D;gene_id=ENSG00000011021;logic_name=ensembl_havana_gene;version=17
1	ensembl_havana	gene	13801445	13840543	.	-	.ID=gene%3AENSG00000162494;Name=LRRC38;biotype=protein_coding;description=leucine rich repeat containing 38 %5BSource%3AHGNC Symbol%3BAcc%3A27005%5D;gene_id=ENSG00000162494;logic_name=ensembl_havana_gene;version=5
1	ensembl_havana	gene	16345370	16360545	.	+	.ID=gene%3AENSG00000186510;Name=CLCNKA;biotype=protein_coding;description=chloride channel%2C voltage-sensitive Ka %5BSource%3AHGNC Symbol%3BAcc%3A2026%5D;gene_id=ENSG00000186510;logic_name=ensembl_havana_gene;version=7
1	ensembl_havana	gene	16370272	16383803	.	+	.ID=gene%3AENSG00000184908;Name=CLCNKB;biotype=protein_coding;description=chloride channel%2C voltage-sensitive Kb %5BSource%3AHGNC Symbol%3BAcc%3A2027%5D;gene_id=ENSG00000184908;logic_name=ensembl_havana_gene;version=13
1	ensembl_havana	gene	25071848	25170815	.	+	.ID=gene%3AENSG00000169504;Name=CLIC4;biotype=protein_coding;description=chloride intracellular channel 4 %5BSource%3AHGNC Symbol%3BAcc%3A13518%5D;gene_id=ENSG00000169504;logic_name=ensembl_havana_gene;version=10
1	ensembl_havana	gene	26517052	26529459	.	+	.ID=gene%3AENSG00000188782;Name=CATSPER4;biotype=protein_coding;description=cation channel%2C sperm associated 4 %5BSource%3AHGNC Symbol%3BAcc%3A23220%5D;gene_id=ENSG00000188782;logic_name=ensembl_havana_gene;version=3
```

END_DOC
 *
 */
@Program(
		name="goutils",
		description="Gene Ontology Utils. Retrieves terms from Gene Ontology",
		keywords={"geneontology","go","gexf"},
		biostars={488538},
		creationDate="20180130",
		modificationDate="20240523",
		jvarkit_amalgamion = true,
		menu="Utilities"
		)
public class GoUtils
	extends Launcher
	{
	private static Logger LOG=Logger.build(GoUtils.class).make();
	private enum Action{dump_table,dump_gexf,goa,gff3};
		
	/** wraps a Go:Term */
	private static class UserTerm
		{
		final GOOntology.Term term;
		final int _hash;
		/** color in GEXF */
		Color vizColor = null;
		/** size in GEXF */
		Double vizSize=null;

		UserTerm(final GOOntology.Term term)
			{
			this.term = term;
			this._hash = term.hashCode();
			}
		@Override
		public int hashCode() {
			return this._hash;
			}
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null) return false;
			return this.term.equals(UserTerm.class.cast(obj).term);
			}
		@Override
		public String toString() {
			return this.term.toString();
			}
		}
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-action","--action"},description=
			"What shoud I do ? default is dump as table. 'goa' only keeps GOA elements in GOA input in GAF format (e.g http://geneontology.org/gene-associations/goa_human.gaf.gz).")
	private Action action = Action.dump_table;
	@Parameter(names="-go",description=GOParser.GO_URL_OPT_DESC)
	private String goURI = GOOntology.GO_OBO_URL;
	@Parameter(names="-goa",description=GOAFileIterator.OPT_DESC)
	private String goaURI = GOAFileIterator.DEFAULT_GOA_URI;
	@Parameter(names={"-g","--gff","--gff3"},description="GFF3 file for action=goa or action=gff3.")
	private String gffPath = null;
	
	@Parameter(names= {"-A","--accession",},description="User Go Terms accession numbers or name.eg GO:0005216 ('ion channel activity') ")
	private Set<String> userAccStrings = new HashSet<>();
	@Parameter(names= {"--exclude-accession","--exclude"},description="User Go Terms to be EXCLUDED accession numbers or name.eg")
	private Set<String> excludeUserAccStrings = new HashSet<>();
	@Parameter(names= {"-af","--accession-file",},
			description="File containing accession numbers. One per line. "
					+ "After the first white space one can define optional attributes for gexf:"
					+ "`color=<COLOR>;size=<SIZE>")
	private Path accessionFile=null;
	@Parameter(names= {"-i","--inverse",},description="inverse the result")
	private boolean inverse=false;
	@Parameter(names= {"--debug",},description="debug",hidden=true)
	private boolean do_debug =false;
    
	private GOOntology mainGoTree=null;
	
	public GoUtils() {
	}
	
	
	
	private void dumpGff3(Gff3Writer out,Gff3Feature feat,final Set<String> geneNames) throws IOException {
		if(feat==null) return;
		if(feat.getType().equals("gene")) {
			if(feat.getAttribute("Name").stream().anyMatch(S->geneNames.contains(S))) {
				out.addFeature(feat);
				for(final Gff3Feature c:feat.getDescendents())  {
					out.addFeature(c);
					}
				}
			}
		else if(feat.hasChildren()) {
			for(Gff3Feature c:feat.getChildren())  {
				dumpGff3(out, c, geneNames);
				}
		}
		}
	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			
			this.mainGoTree = new GOParser().
					setDebug(this.do_debug).
					parseOBO(this.goURI);
			
			final Map<GOOntology.Term,UserTerm> userTerms = new HashMap<>();
			final Map<GOOntology.Term,UserTerm> excludeUserTerms = new HashMap<>();
			for(final String s:this.userAccStrings.stream().flatMap(S->Arrays.stream(S.split("[ ,;]"))).collect(Collectors.toSet()))
				{
				if(StringUtil.isBlank(s)) continue;
				final GOOntology.Term t= this.mainGoTree.getTermByAccessionOrName(s);
				if(t==null)
					{
					LOG.error("cannot find user term \""+s+"\"");
					return -1;
					}
				userTerms.put(t,new UserTerm(t));
				}
			for(final String s:this.excludeUserAccStrings.stream().flatMap(S->Arrays.stream(S.split("[ ,;]"))).collect(Collectors.toSet()))
				{
				if(StringUtil.isBlank(s)) continue;
				final GOOntology.Term t= this.mainGoTree.getTermByAccessionOrName(s);
				if(t==null)
					{
					LOG.error("cannot find user excluded term \""+s+"\"");
					return -1;
					}
				excludeUserTerms.put(t,new UserTerm(t));
				}
			
			
			final Predicate<GOOntology.Term> keepTerm = T->{
				boolean keep=false;
				if(userTerms.isEmpty())
					{
					keep=true;
					}
				else if(userTerms.keySet().
						stream().
						anyMatch(USERTERM->(T.isDescendantOf(USERTERM)))) {
					keep=true;
					}
				if(keep && !excludeUserTerms.isEmpty()) {
					 if(excludeUserTerms.keySet().
								stream().
								anyMatch(USERTERM->(T.isDescendantOf(USERTERM)))) {
						keep = false;
					 	}
					}
				
				if(this.inverse) keep=!keep;
				return keep;
				};
			
			if(this.accessionFile!=null)
				{
				final ColorUtils colorUtils = new ColorUtils();
				try(BufferedReader r=IOUtils.openPathForBufferedReading(this.accessionFile)) {
					String line;
					while((line=r.readLine())!=null) {
						if(line.isEmpty() || line.startsWith("#")) continue;
						int last=0;
						for(last=0;last< line.length();++last) {
							if(Character.isWhitespace(line.charAt(last))) break;
							}
						final String s= line.substring(0, last);
						GOOntology.Term t=  this.mainGoTree.getTermByAccessionOrName(s);
						if(t==null)
								{
								LOG.error("In "+this.accessionFile+" cannot find user term \""+s+"\"");
								return -1;
								}
						final UserTerm ut = new UserTerm(t);
						userTerms.put(t,ut);
						switch(this.action)
							{
						
							case dump_gexf:
									{
									for(final String left: line.substring(last).trim().split("[ \t;]+"))
										{
										if(left.isEmpty())
											{
											// cont
											}
										else if(left.startsWith("color=") && ut.vizColor==null)
											{
											ut.vizColor = colorUtils.parse(left.substring(6));
											}
										else if(left.startsWith("size=") && ut.vizSize==null)
											{
											ut.vizSize = Double.parseDouble(left.substring(5));
											}
										else
											{
											LOG.warning("Ignoring unknown modifier "+left +" in "+line);
											}
										}
									break;
									}
							default:break;
							}
					
						}
					}
				}
			
			
			switch(this.action)
				{
				case dump_gexf:
						{
						final XMLOutputFactory xof=XMLOutputFactory.newFactory();
						XMLStreamWriter w= null;
						FileWriter fw=null;
						
						if(this.outputFile==null)
							{
							w=xof.createXMLStreamWriter(stdout(), "UTF-8");
							}
						else
							{
							w=xof.createXMLStreamWriter((fw=new FileWriter(this.outputFile)));
							}
						final Function<GOOntology.Term,String> term2str=T->T.getAcn().replaceAll("[\\:_#]+", "_");
						w.writeStartDocument("UTF-8", "1.0");
						w.writeStartElement("gexf");
						w.writeAttribute("xmlns", GexfConstants.XMLNS);
						w.writeAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
						w.writeAttribute("xmlns:viz", GexfConstants.XMLNS_VIZ);
						w.writeAttribute("xsi:schemaLocation",GexfConstants.XSI_SCHEMA_LOCATION);
						w.writeAttribute("version", GexfConstants.VERSION);
						w.writeStartElement("meta");
						  w.writeStartElement("creator");
						  w.writeCharacters(getClass().getName()+" by Pierre Lindenbaum");
						  w.writeEndElement();
						 
						  w.writeStartElement("description");
						  w.writeCharacters("Gene Ontology Tree to Gexf :"+getProgramCommandLine());
						  w.writeEndElement();
						
						w.writeEndElement();//meta
						  w.writeStartElement("graph");
						  w.writeAttribute("mode", "static");
						  w.writeAttribute("defaultedgetype", "directed");
						  
						  w.writeStartElement("attributes");
						  w.writeAttribute("class", "edge");
						  w.writeAttribute("mode", "static");
						  w.writeEndElement();//attributes
							
						  w.writeStartElement("attributes");                                                                                     
						  w.writeAttribute("class", "node");
						  w.writeAttribute("mode", "static");
							  
				          w.writeEmptyElement("attribute");
							w.writeAttribute("id", "0");
							w.writeAttribute("title", "description");
							w.writeAttribute("type", "string");
					      
						  w.writeEmptyElement("attribute");
							w.writeAttribute("id", "1");
							w.writeAttribute("title", "accession");
							w.writeAttribute("type", "string");

					      w.writeEmptyElement("attribute");
							w.writeAttribute("id", "2");
							w.writeAttribute("title", "userTerm");
							w.writeAttribute("type", "boolean");
	
					      w.writeEmptyElement("attribute");
							w.writeAttribute("id", "3");
							w.writeAttribute("title", "parentOfUserTerm");
							w.writeAttribute("type", "boolean");
						
					      w.writeEmptyElement("attribute");
							w.writeAttribute("id", "4");
							w.writeAttribute("title", "childOffUserTerm");
							w.writeAttribute("type", "boolean");

					      w.writeEmptyElement("attribute");
							w.writeAttribute("id", "5");
							w.writeAttribute("title", "division");
							w.writeAttribute("type", "boolean");
							
							
						 w.writeEndElement();//attributes
				          
				          
						  
				
						  w.writeStartElement("nodes");
						  w.writeAttribute("count",String.valueOf(this.mainGoTree.size()));
						  for(final GOOntology.Term term:this.mainGoTree.getTerms())
						 	{
							final UserTerm ut = userTerms.get(term);
							w.writeStartElement("node");
							w.writeAttribute("id",term2str.apply(term));
							
							w.writeAttribute("label",term.getName());
							
							w.writeStartElement("attvalues");
							
							  w.writeEmptyElement("attvalue");
								w.writeAttribute("for", "0");
								w.writeAttribute("value",term.getDefinition());
	
							  w.writeEmptyElement("attvalue");
								w.writeAttribute("for", "1");
								w.writeAttribute("value",term.getAcn());
								
							  w.writeEmptyElement("attvalue");
								w.writeAttribute("for", "2");
								w.writeAttribute("value",String.valueOf(ut!=null));
							
							  w.writeEmptyElement("attvalue");
								w.writeAttribute("for", "3");//is parent of any user term
								w.writeAttribute("value",String.valueOf(userTerms.keySet().stream().anyMatch(T->T.isDescendantOf(term))));
							  
							w.writeEmptyElement("attvalue");
								w.writeAttribute("for", "4");//is child of any user term
								w.writeAttribute("value",String.valueOf(userTerms.keySet().stream().anyMatch(T->term.isDescendantOf(T))));

							w.writeEmptyElement("attvalue");
								w.writeAttribute("for", "5");
								w.writeAttribute("value",term.getDivision()==null?".":term.getDivision().name());
							
							w.writeEndElement();//attvalues
							
							double viz_size = 1.0;
							if(ut!=null) {
								if(ut.vizSize!=null)
									{
									viz_size = ut.vizSize;
									}
								if(ut.vizColor!=null)
									{
									// viz:color
									w.writeEmptyElement("viz:color");
										w.writeAttribute("r",String.valueOf(ut.vizColor.getRed()));
										w.writeAttribute("g",String.valueOf(ut.vizColor.getGreen()));
										w.writeAttribute("b",String.valueOf(ut.vizColor.getBlue()));
										w.writeAttribute("a",String.valueOf("1.0"));
									}
								}
							w.writeEmptyElement("viz:size");
								w.writeAttribute("value",String.valueOf(viz_size));
							w.writeEndElement();//node
						 	}
						  w.writeEndElement();//nodes
						
						  
						  w.writeStartElement("edges");
						  w.writeAttribute("count",String.valueOf(this.mainGoTree.getTerms().stream().
								  mapToInt(N->N.getRelations().size()).sum()
								  ));
						  for(final GOOntology.Term term: this.mainGoTree.getTerms())
						 	{
							for(final GOOntology.Relation rel: term.getRelations())
								{
								w.writeStartElement("edge");
								w.writeAttribute("id","E"+term2str.apply(term)+"_"+term2str.apply(rel.getTo()));
								w.writeAttribute("type","directed");
								w.writeAttribute("source",term2str.apply(term));
								w.writeAttribute("target",term2str.apply(rel.getTo()));
								w.writeAttribute("label",rel.getType());
								w.writeAttribute("weight",String.valueOf(1));
								
								final Color vizColor = Color.BLACK;
								
								
								// viz:color
								w.writeEmptyElement("viz:color");
									w.writeAttribute("r",String.valueOf(vizColor.getRed()));
									w.writeAttribute("g",String.valueOf(vizColor.getGreen()));
									w.writeAttribute("b",String.valueOf(vizColor.getBlue()));
									w.writeAttribute("a",String.valueOf("1.0"));
									
								
								w.writeEndElement();
								}
						 	}
						  w.writeEndElement();//edges
						  
						  w.writeEndElement();//graph
				
						
						w.writeEndElement();//gexf
						w.writeEndDocument();
						w.flush();
						if(fw!=null)
							{
							fw.flush();
							CloserUtil.close(fw);
							}
						else
							{
							System.out.flush();
							}
					break;
					}
				case goa:
					{
					if(!args.isEmpty())
						{
						LOG.error("too many arguments");
						return -1;
						}
					final String input;
					if(StringUtil.isBlank(this.goaURI)) {
						input = oneFileOrNull(args);
						}
					else {
						input = this.goaURI;
						}
					final Set<String> acns_set = this.mainGoTree.getTerms().
						stream().
						filter(keepTerm).
						map(T->T.getAcn()).
						collect(Collectors.toSet());
						
					try(BufferedReader br= IOUtils.openURIForBufferedReading(input)) {
						try(GOAFileIterator goain = GOAFileIterator.newInstance(br)) {
							try(PrintWriter out = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
								while(goain.hasNext()) {
									final GOAFileIterator.GafRecord rec = goain.next();
									if(rec.getQualifiers().contains("NOT")) continue;
									if(!acns_set.contains(rec.getGoId())) continue;
									out.println(rec.toString());
									}
								out.flush();
								}
							}
						}
				
					break;
					}
				case gff3:
					{
					if(!args.isEmpty())
						{
						LOG.error("too many arguments");
						return -1;
						}
					if(StringUtil.isBlank(this.goaURI)) {
						LOG.error("undefined GOA-URI");
						return -1;
						}
					final String input;
					if(!StringUtils.isBlank(this.gffPath)) {
						input = oneFileOrNull(args);
						}
					else {
						input = this.gffPath;
						}
					final Set<String> acns_set = this.mainGoTree.getTerms().
						stream().
						filter(keepTerm).
						map(T->T.getAcn()).
						collect(Collectors.toSet());
					final Set<String> geneNames = new HashSet<>();
					try(BufferedReader br= IOUtils.openURIForBufferedReading(this.goaURI)) {
						try(GOAFileIterator goain = GOAFileIterator.newInstance(br)) {
							while(goain.hasNext()) {
								final GOAFileIterator.GafRecord rec = goain.next();
								if(rec.getQualifiers().contains("NOT")) continue;
								if(!acns_set.contains(rec.getGoId())) continue;
								geneNames.add(rec.getObjectSymbol());
								}
							}
						}
					
					final Gff3Codec gff3 = new Gff3Codec(DecodeDepth.DEEP);
					try(InputStream is =(input == null? stdin():IOUtils.openURIForReading(input))) {
						final AsciiLineReader asciiLineReader = AsciiLineReader.from(is);
						final LineIterator lr= new LineIteratorImpl(asciiLineReader);
						try(OutputStream out = super.openFileOrStdoutAsStream(this.outputFile)) {
							Gff3Writer gw = new Gff3Writer(out);
							while(!gff3.isDone(lr)) {
								dumpGff3(gw,gff3.decode(lr),geneNames);
								}
							out.flush();
							}
						gff3.close(lr);
						asciiLineReader.close();
						}
					break;
					}
				case dump_table://through
				default:
					{
					if(!args.isEmpty()) {
						LOG.error("too many arguments");
						return -1;
						}
					try(PrintWriter out = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
						out.println("#ACN\tNAME\tDEFINITION\tDIVISION");
						for(final GOOntology.Term t:this.mainGoTree.getTerms())
							{
							if(keepTerm.test(t))
								{
								out.print(t.getAcn());
								out.print('\t');
								out.print(t.getName());
								out.print('\t');
								out.print(t.getDefinition());
								out.print('\t');
								out.print(t.getDivision()==null?".":t.getDivision().name());
								out.println();
								}
							}
						out.flush();
						}
					break;
					}
				}
				return 0;
				}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}

	public static void main(final String[] args)
		{
		new GoUtils().instanceMainWithExit(args);
		}
	}
