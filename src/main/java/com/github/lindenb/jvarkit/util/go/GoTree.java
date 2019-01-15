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

*/
package com.github.lindenb.jvarkit.util.go;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.RDF;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;

/**
 * Gene Ontology tree
 * @author lindenb
 *
 */
public class GoTree implements Iterable<GoTree.Term>
	{
	private static Logger LOG=Logger.build(GoTree.class).make(); 

	private static final String NS="http://www.geneontology.org/dtds/go.dtd#";
	public static final String GO_RDF_URL="http://archive.geneontology.org/latest-termdb/go_daily-termdb.rdf-xml.gz";
	public static final String GO_URL_OPT_DESC="Gene ontology URI. Formatted as RDF+XML. Can be gzipped.";
	
	public static enum RelType
		{
		is_a,negatively_regulates,part_of,positively_regulates,regulates
		}
	private static final Set<RelType> REL_TYPE_SET = Arrays.asList(RelType.values()).stream().collect(Collectors.toSet());

	public static enum Division
		{
		cellular_component("GO:0005575"),biological_process("GO:0008150"),molecular_function("GO:0003674");
		private final String acn;
		Division(final String acn) { this.acn=acn;}
		/** get accession number */
		public String getAcn() { return this.acn;}
		}
	private static final Set<Division> DIVISION_SET = Arrays.asList(Division.values()).stream().collect(Collectors.toSet());

	
	public static interface Relation
		{
		public RelType getType();
		public Term getTo();
		}
	
	public static class RelTypeConverter implements IStringConverter<Set<RelType>>
		{
		@Override
		public Set<RelType> convert(final String s) {
			if(StringUtil.isBlank(s))
				{
				return REL_TYPE_SET;
				}
			final Set<RelType> rt = new HashSet<>();
			for(final String token:s.split("[,; |]"))
				{
				if(StringUtil.isBlank(token)) continue;
				final Optional<RelType> f = REL_TYPE_SET.stream().
						filter(T->T.name().equalsIgnoreCase(token)).
						findFirst();
				if(!f.isPresent())
					{
					throw new JvarkitException.UserError(
						"undefined GO:rel "+token+" available are: " +
						REL_TYPE_SET
						);
					}
				
				rt.add(f.get());
				}
			if(rt.isEmpty()) return REL_TYPE_SET;
			return rt;
			}
		}
	
	
	public static class DivisionConverter implements IStringConverter<Set<Division>>
		{
		@Override
		public Set<Division> convert(final String s) {
			if(StringUtil.isBlank(s))
				{
				return DIVISION_SET;
				}
			final Set<Division> ds = new HashSet<>();
			for(final String token:s.split("[,; |]"))
				{
				if(StringUtil.isBlank(token)) continue;
				final Optional<Division> f = DIVISION_SET.stream().
						filter(T->T.name().equalsIgnoreCase(token)).
						findFirst();
				if(!f.isPresent())
					{
					throw new JvarkitException.UserError(
						"undefined division "+token+" available are: " +
								DIVISION_SET
						);
					}
				
				ds.add(f.get());
				}
			if(ds.isEmpty()) return DIVISION_SET;
			return ds;
			}
		}

	
	public static class ReadingGo
		{
		@Parameter(names={"-go","--go","--gene-ontology"},description=GoTree.GO_URL_OPT_DESC)
		public String goUri = GoTree.GO_RDF_URL;
		@Parameter(names={"-go-relations","--go-relations"},description="limit the gene ontology tree to those relationships. empty: all possible relationships. ",converter=RelTypeConverter.class)
		public Set<RelType> relTypes = REL_TYPE_SET;
		@Parameter(names={"-go-divisions","--go-divisions"},description="limit the gene ontology tree to those divisions. empty: all possible divisions. ",converter=DivisionConverter.class)
		public Set<Division> divisions = DIVISION_SET;

		public Parser createParser()
			{
			return new Parser().
					setDivisions(this.divisions).
					setRelations(this.relTypes);
			}
		}
	
	private static class RelationImpl implements Relation
		{
		final TermImpl termFrom;
		final RelType reltype;
		final TermImpl termTo;
		final private int _hash;
		RelationImpl(final TermImpl termFrom,final RelType reltype,final TermImpl termTo)
			{
			this.termFrom = termFrom;
			this.reltype = reltype;
			this.termTo = termTo;
			
			final int prime = 31;
			int result = 1;
			result = prime * result + termFrom.hashCode();
			result = prime * result + reltype.hashCode();
			result = prime * result + termTo.hashCode();
			this._hash = result;
			}
		@Override
		public Term getTo() {
			return termTo;
			}
		@Override
		public RelType getType() {
			return reltype;
			}
		@Override
		public int hashCode() {
			return this._hash;
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			final RelationImpl other = (RelationImpl) obj;
			return (reltype == other.reltype &&
					termFrom.equals(other.termFrom) &&
					termTo.equals(other.termTo)
					);
			}
		@Override
		public String toString() {
			return this.termFrom.getAcn()+"["+this.termFrom.getName()+"] -- "+
					this.reltype.name()+ " --> "+
					this.termTo.getAcn()+"["+this.termTo.getName()
					;
			}
		}
	
	public static interface Term
		{
		/** get accession number */
		public String getAcn();
		/** return name like 'maltose catabolic process' GO:0000025 */
		public String getName();
		/** return true if getRelations() is not empty */
		public boolean hasRelations();
		
		public Set<String> getSynonyms();
		public String getDefinition();
		public Set<Relation> getRelations();
		public boolean isDescendantOf(Term other);
		public List<DbXRef> getDbXRefs();
		/** get Min depth, root/all is '0' */
		public int getMinDepth();
		}
	
	public static interface DbXRef
		{
		public String getDatabaseSymbol();
		public String getReference();
		}
	private static class DbXRefImpl implements DbXRef
		{
		private String symbol=null;
		private String reference=null;
		public String getDatabaseSymbol()
			{
			return symbol;
			}
		public String getReference()
			{
			return reference;
			}
		}
	
	
	private GoTree()
		{
		}
	
	public int size()
		{
		return acn2term.size();
		}
	
	public Collection<Term> getTerms() 
		{
		return Collections.unmodifiableCollection(this.acn2term.values());
		}
	
	@Override
	public Iterator<Term> iterator() {
		return getTerms().iterator();
		}
	
	
	private final HashMap<String, TermImpl> acn2term=new HashMap<String, TermImpl>();
	
	private class TermImpl implements Term
		{
		String accession;
		String name;
		String definition = null;
		Set<String> synonyms = null;
		List<DbXRef> dbxrefs = null;
		//final Set<TermImpl> parents=new HashSet<>();
		//final Set<TermImpl> children=new HashSet<>();
		final Set<RelationImpl> relations = new HashSet<>();
		
		TermImpl(final String uri) {
			int hash=uri.indexOf('#');
			if(hash!=-1)
				{
				this.accession=uri.substring(hash+1);
				}
			else
				{
				this.accession = uri;
				}
			this.name= this.accession;
			}
		@Override
		public Set<String> getSynonyms()
			{
			return this.synonyms==null?
					Collections.emptySet():
					Collections.unmodifiableSet(this.synonyms)
					;
			}
		
		@Override
		public String getAcn() {
			return accession;
			}
		@Override
		public String getName() {
			return name;
			}
		@Override
		public String getDefinition() {
			return StringUtil.isBlank(definition)?getName():this.definition;
			}
		@Override
		public List<DbXRef> getDbXRefs()
			{
			return this.dbxrefs==null?
					Collections.emptyList():
					Collections.unmodifiableList(this.dbxrefs)
					;
			}
		@Override
		public boolean hasRelations()
			{
			return !this.relations.isEmpty();
			}
		
		@Override
		public Set<Relation> getRelations() {
			return Collections.unmodifiableSet(this.relations);
			}
		
		@Override
		public boolean isDescendantOf(final Term parentNode)
			{
			if(parentNode==this || getAcn().equals(parentNode.getAcn())) return true;
			for(final RelationImpl rel:this.relations)
				{
				if(rel.getTo().isDescendantOf(parentNode)) return true;
				}
			return false;
			}
		
		@Override
		public int hashCode()
			{
			return this.accession.hashCode();
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj) return true;
			if (obj == null) return false;
			if (getClass() != obj.getClass()) return false;
			final TermImpl other = (TermImpl) obj;
			return this.accession.equals(other.accession);
			}
		
		public int getMinDepth()
			{
			if(!hasRelations()) return 0;
			int d=-1;
			for(final RelationImpl rel:this.relations)
				{
				int d2 = 1 + rel.getTo().getMinDepth();
				if(d==-1 || d>d2)
					{
					d=d2;
					}
				}
			return d;
			}
		
		@Override
		public String toString() {
			return accession;
			}
		}
	
	public Term getTermByAccession(final String s)
		{
		return this.acn2term.get(s);
		}
	
	/** search term by name or accession or synonyms, ignoring case */
	public Term getTermByName(final String s)
		{
		final Term t0 = getTermByAccession(s);
		if(t0!=null) return t0;
		for(final TermImpl t:this.acn2term.values())
			{
			if(t.accession.equalsIgnoreCase(s)) return t;
			if(t.name!=null && t.name.equalsIgnoreCase(s)) return t;
			if(t.synonyms!=null)
				{
				for(final String syn:t.synonyms)
					{
					if(syn.equalsIgnoreCase(s)) return t;
					}
				}
			}
		return null;
		}
		
	@Deprecated
	public static GoTree parse(Reader xmlIn) throws IOException,XMLStreamException
		{
		return new Parser().parse(xmlIn);
		}
	@Deprecated
	public static GoTree parse(File file) throws IOException,XMLStreamException
		{
		return new Parser().parse(file);
		}
	@Deprecated
	public static GoTree parse(String uri) throws IOException,XMLStreamException
		{
		return new Parser().parse(uri);
		}
	
	public static class Parser
		{
		private static final QName rdfAbout=new QName(RDF.NS,"about",RDF.pfx);
		private static final QName rdfRsrc=new QName(RDF.NS,"resource",RDF.pfx);
		private static final QName parseType=new QName(RDF.NS,"parseType",RDF.pfx);
		
		private final Map<String,TermImpl> uri2term = new HashMap<>();
		private final List<IsA> isAList=new ArrayList<>();
		private boolean ignore_synonyms=false;
		private boolean ignore_definitions=false;
		private boolean ignore_dbxref=false;
		private boolean debug=false;
		private Set<RelType> userRelTypes = GoTree.REL_TYPE_SET;
		private Set<Division> userDivisions = GoTree.DIVISION_SET;
		
		
		/** set accepted relations */
		public Parser setRelations(final Set<RelType> userRelTypes) {
			this.userRelTypes = (userRelTypes==null? GoTree.REL_TYPE_SET:new HashSet<>(userRelTypes));
			return this;
			}
		/** set accepted divisions */
		public Parser setDivisions(final Set<Division> userRelTypes) {
			this.userDivisions = (userRelTypes==null? GoTree.DIVISION_SET:new HashSet<>(userRelTypes));
			return this;
			}
		
		public Parser setIgnoreDbXRef(boolean ignore_dbxref) {
			this.ignore_dbxref = ignore_dbxref;
			return this;
			}
		
		public Parser setIgnoreDefinitions(boolean ignore_definitions) {
			this.ignore_definitions = ignore_definitions;
			return this;
		}
		
		public Parser setIgnoreSynonyms(boolean ignore_synonyms) {
			this.ignore_synonyms = ignore_synonyms;
			return this;
			}
		
		public Parser setDebug(boolean d) {
			this.debug = d;
			return this;
			}
		
		private static class IsA
			{
			final TermImpl term;
			final String parentUri;
			final RelType relType;
			IsA(final TermImpl term,final String parentUri,final RelType relType)
				{
				this.term = term;
				this.parentUri = parentUri;
				this.relType=relType;
				}
			}
		
		private DbXRef parseDbXRefRDF(
				final StartElement root,
				final XMLEventReader r) throws IOException,XMLStreamException
			{
			final Attribute parseTypeAtt=root.getAttributeByName(parseType);
			if(parseTypeAtt==null)
				{
				throw new XMLStreamException("no rdf:parseType",root.getLocation());
				}
			if(!parseTypeAtt.getValue().equals("Resource")) {
				throw new XMLStreamException("no rdf:parseType=Resource",root.getLocation());
				}
			final DbXRefImpl xref = new DbXRefImpl();
			while(r.hasNext())
				{
				XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					final StartElement E=evt.asStartElement();
					final QName qN=E.getName();
					if( NS.equals(qN.getNamespaceURI()))
						{
						final String localName = qN.getLocalPart();
						if(localName.equals("database_symbol"))
							{
							xref.symbol=r.getElementText();
							}
						else if(localName.equals("reference"))
							{
							xref.reference=r.getElementText();
							}
						}
					
					}
				else if(evt.isEndElement())
					{
					final EndElement E=evt.asEndElement();
					final QName qN=E.getName();
					if(qN.getLocalPart().equals("dbxref") && NS.equals(qN.getNamespaceURI()))
						{
						break;
						}
					}
				}
			if(StringUtil.isBlank(xref.symbol)) return null;
			if(StringUtil.isBlank(xref.reference)) return null;
			return xref;
			}
		
		
		private RelType findRelTypeByName(final String s)
			{
			for(RelType rt:this.userRelTypes)
				{
				if(s.equals(rt.name())) return rt;
				}
			return null;
			}
		
		private void parseRDFTerm(
				final GoTree tree,
				final StartElement root,
				final XMLEventReader r) throws IOException,XMLStreamException
			{
			final Attribute aboutAtt=root.getAttributeByName(rdfAbout);
			if(aboutAtt==null)
				{
				throw new XMLStreamException("no rdf:about",root.getLocation());
				}
			final String termUri=aboutAtt.getValue();
			if(this.debug) LOG.debug("found term "+termUri);
			/* term exists  ? */
			if(this.uri2term.containsKey(termUri))
				{
				throw new XMLStreamException("Term URI defined twice:"+termUri,root.getLocation());
				}
			boolean obsolete_flag=false;
			final TermImpl term = tree.new TermImpl(termUri);
				
			while(r.hasNext())
				{
				XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					final StartElement E=evt.asStartElement();
					final QName qN=E.getName();
					if( NS.equals(qN.getNamespaceURI()))
						{
						final String localName = qN.getLocalPart();
						RelType reltype=null;
						if(localName.equals("accession"))
							{
							term.accession=r.getElementText();
							}
						else if(localName.equals("name"))
							{
							term.name=r.getElementText();
							}
						else if(localName.equals("synonym") && !this.ignore_synonyms)
							{
							if(term.synonyms==null) term.synonyms=new HashSet<>();
							term.synonyms.add(r.getElementText());
							}
						else if(localName.equals("definition") && !this.ignore_definitions)
							{
							term.definition=r.getElementText();
							}
						else if((reltype = findRelTypeByName(localName))!=null)
							{
							final Attribute rsrc=E.getAttributeByName(rdfRsrc);
							if(rsrc==null) throw new XMLStreamException("att missing "+rdfRsrc+" for "+aboutAtt.getValue(),evt.getLocation());
							final String parentUri=rsrc.getValue();
							if(parentUri.startsWith("http://www.geneontology.org/go#obsolete")) {
								obsolete_flag=true;
								}
							if(!obsolete_flag)
								{
								this.isAList.add(new IsA(term, parentUri,reltype));
								}
							}
						else if(localName.equals("dbxref"))
							{
							final DbXRef xref = parseDbXRefRDF(E,r);
							if(!this.ignore_dbxref && xref!=null && !obsolete_flag)
								{
								if(term.dbxrefs==null) term.dbxrefs=new ArrayList<>();
								term.dbxrefs.add(xref);
								}
							}
						}
					
					}
				else if(evt.isEndElement())
					{
					final EndElement E=evt.asEndElement();
					final QName qN=E.getName();
					if(qN.getLocalPart().equals("term") && 
					  NS.equals(qN.getNamespaceURI()))
						{
						/* found this, no name, linked nowhere... */
						if(StringUtil.isBlank(term.name))
							{
							obsolete_flag = true;
							}
						
						if(!obsolete_flag) 
							{
							this.uri2term.put(termUri, term);
							}
						break;
						}
					}
				}
		
			}	
		
		public GoTree parse(final Reader xmlIn) throws IOException,XMLStreamException
			{
			final GoTree tree=new GoTree();
			XMLInputFactory fact=XMLInputFactory.newFactory();
			fact.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
			XMLEventReader r=fact.createXMLEventReader(xmlIn);
			while(r.hasNext())
				{
				final XMLEvent evt=r.nextEvent();
				if(evt.isStartElement())
					{
					final StartElement E=evt.asStartElement();
					final QName qN=E.getName();
					if(qN.getLocalPart().equals("term") && 
						NS.equals(qN.getNamespaceURI()))
						{
						this.parseRDFTerm(tree,E,r);
						}
					}
				}
			r.close();
			for(final TermImpl t:this.uri2term.values()) {
				tree.acn2term.put(t.getAcn(), t);
				}
			
			for(final IsA isa:this.isAList)
				{
				final TermImpl parentTerm = this.uri2term.get(isa.parentUri);
				if(parentTerm==null) {
					LOG.warning("Cannot find uri :"+isa.parentUri);
					continue;
					}
				if(isa.term.equals(parentTerm))
					{
					LOG.warning("self/self relation ??"+isa.parentUri);
					continue;
					}
				isa.term.relations.add(
						new RelationImpl(isa.term, isa.relType,parentTerm)
						);
				
				}
			
			if(!this.userDivisions.equals(DIVISION_SET))
				{
				final Set<String> toRemove=new HashSet<>();
				for(final Term term:tree.acn2term.values())
					{
					boolean keep=false;
					for(final Division div:this.userDivisions)
						{
						final Term root= tree.getTermByAccession(div.getAcn());
						if(root==null) throw new JvarkitException.UserError("cannot find "+div+" in go tree !");
						if(root.equals(term) ||root.isDescendantOf(term) /* 'all'*/ || term.isDescendantOf(root))
							{
							keep=true;
							break;
							}
						}
					if(!keep) toRemove.add(term.getAcn());
					}
				for(final String acn:toRemove) tree.acn2term.remove(acn);
				}
			
			if(this.debug)
				{
				LOG.debug("tree size: "+tree.acn2term.size());
				}
			this.uri2term.clear();
			this.isAList.clear();
			return tree;
			}
		
		
		public GoTree parse(final File file) throws IOException,XMLStreamException
			{
			Reader r=null;
			try
				{
				r=IOUtils.openFileForReader(file);
				return parse(r);
				}
			finally
				{
				CloserUtil.close(r);
				}
			}
		public GoTree parse(final String uri) throws IOException,XMLStreamException
			{
			Reader r=null;
			try
				{
				r=new InputStreamReader(IOUtils.openURIForReading(uri),"UTF-8");
				return parse(r);
				}
			finally
				{
				CloserUtil.close(r);
				}
			}
		}
	
	}
