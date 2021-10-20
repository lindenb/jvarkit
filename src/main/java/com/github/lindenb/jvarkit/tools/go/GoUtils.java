/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.io.PrintWriter;
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
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.gexf.GexfConstants;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.go.GoTree;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.swing.ColorUtils;
import com.github.lindenb.jvarkit.goa.GOAFileIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

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

Use GO annotation to retrieve genes associated to GO:0005216 'ion channel activity' 

```
join -t $'\t' -1 1 -2 2 \
	<(java -jar dist/goutils.jar -A 'GO:0005216' | cut -f 1 | sort | uniq) \
	<(wget -q -O - "http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD" | gunzip -c | grep -v '^!' | cut -f3,5 | uniq | LC_ALL=C sort -t $'\t' -k2,2) |\
sort -t $'\t' -k2,2 |\
grep SCN5A -A 10 -B 10
```

```
(...)
GO:0086006	SCN2B
GO:0005244	SCN3A
GO:0005248	SCN3A
GO:0005248	SCN3A
GO:0005248	SCN3B
GO:0086006	SCN3B
GO:0005248	SCN4A
GO:0005248	SCN4A
GO:0005248	SCN4B
GO:0086006	SCN4B
GO:0005244	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0005248	SCN5A
GO:0086006	SCN5A
GO:0086060	SCN5A
GO:0086061	SCN5A
GO:0086062	SCN5A
GO:0086063	SCN5A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN7A
GO:0005248	SCN8A
GO:0005248	SCN8A
GO:0005248	SCN9A
GO:0005248	SCN9A
GO:0005248	SCN9A
GO:0005272	SCNN1A
(...)
```

```
$ wget -q -O - "http://geneontology.org/gene-associations/goa_human.gaf.gz" | gunzip -c | java -jar dist/goutils.jar -go go.rdf.gz --action goa -A 'ion transmembrane transport'  | head
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

END_DOC
 *
 */
@Program(
		name="goutils",
		description="Gene Ontology Utils. Retrieves terms from Gene Ontology",
		keywords={"geneontology","go","gexf"},
		biostars={488538},
		creationDate="20180130",
		modificationDate="20211020"
		)
public class GoUtils
	extends Launcher
	{
	private static Logger LOG=Logger.build(GoUtils.class).make();
	private enum Action{dump_table,dump_gexf,goa};
	
	/** wraps a Go:Term */
	private static class UserTerm
		{
		final GoTree.Term term;
		final int _hash;
		/** color in GEXF */
		Color vizColor = null;
		/** size in GEXF */
		Double vizSize=null;

		UserTerm(final GoTree.Term term)
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
	@ParametersDelegate
	private GoTree.ReadingGo readingGo = new GoTree.ReadingGo();
	@Parameter(names= {"-A","--accession",},description="User Go Terms accession numbers or name")
	private Set<String> userAccStrings = new HashSet<>();
	@Parameter(names= {"-af","--accession-file",},
			description="File containing accession numbers. One per line. "
					+ "After the first white space one can define optional attributes for gexf:"
					+ "`color=<COLOR>;size=<SIZE>")
	private File accessionFile=null;
	@Parameter(names= {"-i","--inverse",},description="inverse the result")
	private boolean inverse=false;
	@Parameter(names= {"--debug",},description="debug",hidden=true)
	private boolean do_debug =false;
    
	private GoTree mainGoTree=null;
	
	public GoUtils() {
	}
	
	private GoTree.Term findTerm(final String s)
		{
		GoTree.Term t= this.mainGoTree.getTermByAccession(s);
		if(t==null)
			{
			t= this.mainGoTree.getTermByName(s);
			}
		return t;
		}
	
	@Override
	public int doWork(final List<String> args) {			
		try
			{
			
			this.mainGoTree = this.readingGo.createParser().
					setDebug(this.do_debug).
					parse(this.readingGo.goUri);
			
			final Map<GoTree.Term,UserTerm> userTerms = new HashMap<>();
			for(final String s:this.userAccStrings)
				{
				if(StringUtil.isBlank(s)) continue;
				final GoTree.Term t= this.findTerm(s);
				if(t==null)
					{
					LOG.error("cannot find user term \""+s+"\"");
					return -1;
					}
				userTerms.put(t,new UserTerm(t));
				}
			final Predicate<GoTree.Term> keepTerm = T->{
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
				if(this.inverse) keep=!keep;
				return keep;
				};
			
			if(this.accessionFile!=null)
				{
				final ColorUtils colorUtils = new ColorUtils();
				try(BufferedReader r=IOUtils.openFileForBufferedReading(this.accessionFile)) {
					String line;
					while((line=r.readLine())!=null) {
						if(line.isEmpty() || line.startsWith("#")) continue;
						int last=0;
						for(last=0;last< line.length();++last) {
							if(Character.isWhitespace(line.charAt(last))) break;
							}
						final String s= line.substring(0, last);
						GoTree.Term t=  this.findTerm(s);
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
						final Function<GoTree.Term,String> term2str=T->T.getAcn().replaceAll("[\\:_#]+", "_");
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
						  for(final GoTree.Term term:this.mainGoTree.getTerms())
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
						  for(final GoTree.Term term: this.mainGoTree.getTerms())
						 	{
							for(final GoTree.Relation rel: term.getRelations())
								{
								w.writeStartElement("edge");
								w.writeAttribute("id","E"+term2str.apply(term)+"_"+term2str.apply(rel.getTo()));
								w.writeAttribute("type","directed");
								w.writeAttribute("source",term2str.apply(term));
								w.writeAttribute("target",term2str.apply(rel.getTo()));
								w.writeAttribute("label",rel.getType().getName());
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
					final Set<String> acns_set = this.mainGoTree.getTerms().
						stream().
						filter(keepTerm).
						map(T->T.getAcn()).
						collect(Collectors.toSet());
						
					final String input = oneFileOrNull(args);
					try(BufferedReader br= super.openBufferedReader(input)) {
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
				case dump_table://through
				default:
					{
					if(!args.isEmpty())
						{
						LOG.error("too many arguments");
						return -1;
						}
					try(PrintWriter out = super.openFileOrStdoutAsPrintWriter(this.outputFile)) {
						out.println("#ACN\tNAME\tDEFINITION\tDIVISION");
						for(final GoTree.Term t:this.mainGoTree.getTerms())
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
