/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.samplesrdf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalInt;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import javax.xml.namespace.QName;
import javax.xml.stream.EventFilter;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.XMLEvent;

import org.apache.jena.rdf.model.Literal;
import org.apache.jena.rdf.model.Model;
import org.apache.jena.rdf.model.ModelFactory;
import org.apache.jena.rdf.model.Property;
import org.apache.jena.rdf.model.RDFNode;
import org.apache.jena.rdf.model.Resource;
import org.apache.jena.rdf.model.ResourceFactory;
import org.apache.jena.rdfxml.xmlinput.StAX2Model;
import org.apache.jena.vocabulary.OWL;
import org.apache.jena.vocabulary.RDF;
import org.apache.jena.vocabulary.RDFS;
import org.xml.sax.SAXException;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jena.JenaUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.stream.CollectorsUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.StringUtil;
/**
BEGIN_DOC

# Motivation

A tool formatting samples definition.

# Input Format

Input is a Recfile ( https://en.wikipedia.org/wiki/Recfiles )  in the form:

```
key : value
key : value
key : value
key : value

key : value
key : value
key : value
key : value

key : value
(...)
```

the key/value delimiter is `:` or `=`.
multiple records are separated by one or more spaces.
lines starting with "#" are ignored.
keys are case insensible.

**type:** : each record requires the key `type`. The type is either Sample or Group. A Sample defines one Sample, A group is a set of samples.

**id:** : each record requires the key `id`. The type is either Sample or Group. It defines a unique identifier in the database. Id are converted to uppercase

**label:** : a title for this record

**comment:** or '**description:** : a description for the record

** hpo: ** : or '**hp:** is a identifier in the human phenotype ontology database. If the `hpo:` is defined in a `type:Group`, the phenotype will be propaged to all the samples in that group.
The following syntaxes are allowed:
```
hp: HP:0011712 everything after the first word is ignored
hp: 0011712 words
hpo: 0011712 words
```


** doid:** is a identifier in the human disease ontology database. If the `doid:` is defined in a `type:Group`, the phenotype will be propaged to all the samples in that group.
The following syntaxes are allowed:
```
doid: DOID:0011712 everything after the first word is ignored
doid: 0011712 words
doid: 0011712 words
```



## Samples only:
**family:** a family name for this a `Sample`

**father:** the id of the father 

**sex:** the sex of the sample 'male' or 'female'

**birth:** the  year of birth

**alias:** another name for this sample


## group:

**group: id ** or **extends:id** : the current group extends another group defined by 'id'


**pop** or **population** is an identifier in the SNOMED population ontology. https://bioportal.bioontology.org/ontologies/SNOMED-Ethnic-Grou . Value is either the label or the SNOMEDID in the OWL ontology
```
type: sample
id: John
pop: Danes
pop: S-61240
```


 
# Output Format

 * (prefix)_errors.txt : errors found
 * (prefix)_model.rdf : output formatted as RDF/XML
 * (prefix)_doid.txt: output each item in Disease Ontology. Each line is a term in DO listing each sample with this disease.
 * (prefix)_pedigree.txt : output samples formatted as a pedigree file. The links / sex are not checked
 * (prefix)_samples.txt : output samples with their diseases, phenotypes...


# Example

```
$ cat input.recfile 

type:sample
id: azd
father: x1
sex: male
doid: DOID:0050451 Hello
hpo: HP:0011712 World
pop: S-61020

type:group
id:g1
sample: azd
sample: x1


 java -jar dist/jvarkit.jar samplesrdf --population  pop.owl --doid doid.owl --hpo hp.owl -o JETER input.recfile
 
 $ more JETER_doid.txt 
http://purl.obolibrary.org/obo/DOID_0050451	DOID:0050451	Brugada syndrome	AZD
http://purl.obolibrary.org/obo/DOID_0110221	DOID:0110221	Brugada syndrome 4	AZD
http://purl.obolibrary.org/obo/DOID_0110222	DOID:0110222	Brugada syndrome 5	AZD
http://purl.obolibrary.org/obo/DOID_0110220	DOID:0110220	Brugada syndrome 3	AZD
http://purl.obolibrary.org/obo/DOID_0110225	DOID:0110225	Brugada syndrome 8	AZD
http://purl.obolibrary.org/obo/DOID_0110226	DOID:0110226	Brugada syndrome 9	AZD
http://purl.obolibrary.org/obo/DOID_0110223	DOID:0110223	Brugada syndrome 6	AZD
http://purl.obolibrary.org/obo/DOID_0110224	DOID:0110224	Brugada syndrome 7	AZD
http://purl.obolibrary.org/obo/DOID_0110219	DOID:0110219	Brugada syndrome 2	AZD
http://purl.obolibrary.org/obo/DOID_0110218	DOID:0110218	Brugada syndrome 1	AZD


$ more JETER_hpo.txt 
http://purl.obolibrary.org/obo/HP_0011712	HP:0011712	Right bundle branch block	AZD

$ more JETER_groups.txt 
G1	G1		2	1	0	AZD;X1


$ more JETER_samples.txt
X1	X1		.	.	.	.	.	.	.	.	G1
AZD	AZD		.	.	X1	.	male	Right bundle branch block	Brugada syndrome 8; Brugada syndrome 7; Brugada syndrome 6; Brugada syndrome 5; Brugada syndrome 1; B
rugada syndrome; Brugada syndrome 9; Brugada syndrome 2; Brugada syndrome 4; Brugada syndrome 3	Armenians	G1

$ more JETER_pedigree.txt 
X1	X1	0	0	0
AZD	AZD	X1	0	1

$ more JETER_errors.txt 
<https://umr1087.univ-nantes.fr/db/X1> rdf:type <https://umr1087.univ-nantes.fr/Sample> was used but not declared
X1 declared as father of AZD. But sex is not male

$ more JETER_model.rdf 
<rdf:RDF
    xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
    xmlns:u="https://umr1087.univ-nantes.fr/" > 
  <rdf:Description rdf:about="https://umr1087.univ-nantes.fr/db/G1">
    <u:contains-sample rdf:resource="https://umr1087.univ-nantes.fr/db/AZD"/>
    <u:id>G1</u:id>
    <rdf:type rdf:resource="https://umr1087.univ-nantes.fr/Group"/>
    <u:contains-sample rdf:resource="https://umr1087.univ-nantes.fr/db/X1"/>
  </rdf:Description>
  <rdf:Description rdf:about="https://umr1087.univ-nantes.fr/db/X1">
    <u:id>X1</u:id>
    <rdf:type rdf:resource="https://umr1087.univ-nantes.fr/Sample"/>
  </rdf:Description>
  <rdf:Description rdf:about="https://umr1087.univ-nantes.fr/db/AZD">
    <u:hpo rdf:resource="http://purl.obolibrary.org/obo/HP_0011712"/>
    <u:doid rdf:resource="http://purl.obolibrary.org/obo/DOID_0050451"/>
    <u:sex rdf:resource="https://umr1087.univ-nantes.fr/Male"/>
    <u:father rdf:resource="https://umr1087.univ-nantes.fr/db/X1"/>
    <u:id>AZD</u:id>
    <rdf:type rdf:resource="https://umr1087.univ-nantes.fr/Sample"/>
    <u:population rdf:resource="http://purl.bioontology.org/ontology/SNOMED-Ethnic-Group#113169009"/>
  </rdf:Description>
</rdf:RDF>
```

END_DOC
 */
@Program(name="samplesrdf",
description="Digests a  database of samples from a set of recfiles",
keywords={"samples","database","rdf","ontology"},
creationDate="20230201",
modificationDate="20230202",
jvarkit_amalgamion = true,
menu="Utilities"
)
public class SamplesRDF extends Launcher {
	private static final Logger LOG = Logger.build(SamplesRDF.class).make();

	private static final String NS = "https://umr1087.univ-nantes.fr/";
	private static final String FOAF_NS = "http://xmlns.com/foaf/0.1/";
	private static final String OBOINOWL = "http://www.geneontology.org/formats/oboInOwl#";
	private static final String SNOMED_NS="http://purl.bioontology.org/ontology/SNOMED-Ethnic-Group#";
	private static final Resource TYPE_SAMPLE = ResourceFactory.createResource(FOAF_NS+"Person");
	private static final Resource TYPE_GROUP = ResourceFactory.createResource(FOAF_NS+"Group");
	/* flag the Sample/Group that have been explicity created */
	private static final Resource TYPE_DECLARED = ResourceFactory.createResource(NS+"_Validated");
	private static final Resource IS_MALE = ResourceFactory.createResource("http://purl.bioontology.org/ontology/SNOMEDCT/248153007");
	private static final Resource IS_FEMALE = ResourceFactory.createResource("http://purl.bioontology.org/ontology/SNOMEDCT/248152002");
	private static final Property PROP_ID = ResourceFactory.createProperty(NS,"id");
	private static final Property PROP_HPO = ResourceFactory.createProperty(NS,"hpo");
	private static final Property PROP_POPULATION = ResourceFactory.createProperty(NS,"population");
	private static final Property PROP_DOID = ResourceFactory.createProperty(NS,"doid");
	private static final Property PROP_alias = ResourceFactory.createProperty(NS,"alias");
	private static final Property PROP_father = ResourceFactory.createProperty(NS,"father");
	private static final Property PROP_mother = ResourceFactory.createProperty(NS,"mother");
	private static final Property PROP_family = ResourceFactory.createProperty(FOAF_NS,"family_name");
	private static final Property PROP_sex = ResourceFactory.createProperty(FOAF_NS,"gender");
	private static final Property PROP_birthYear = ResourceFactory.createProperty(NS,"birthYear");
	private static final Property OBOINOWBL_id = ResourceFactory.createProperty(OBOINOWL,"id");
	private static final Property SNOMEDID = ResourceFactory.createProperty(SNOMED_NS,"SNOMEDID");
	private static final Property PROP_contains_sample = ResourceFactory.createProperty(FOAF_NS,"member");
	private static final Property PROP_group = ResourceFactory.createProperty(NS,"group");

	
	@Parameter(names={"-hpo","--hpo"},description="OWL-formatted HPO ontology https://github.com/obophenotype/human-phenotype-ontology/blob/master/hp.owl")
	private Path hpoFilename = null;
	@Parameter(names={"-doid","--doid"},description="OWL-formatted Disease ontology. https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.owl ")
	private  Path doidFilename = null;
	@Parameter(names={"-population","--population","--pop","--snomed"},description="OWL-formatted SNOWMED-ETHNIC-GROUP https://bioportal.bioontology.org/ontologies/SNOMED-Ethnic-Group")
	private  Path snomedEthnicFilename = null;

	@Parameter(names={"-o","--out"},description="Base filename for output files",required = true)
	protected String outputBase = null;


	private final Model model;
	
	
	
	public SamplesRDF() {
		this.model = ModelFactory.createDefaultModel();
		this.model.setNsPrefix("u",NS);
		this.model.setNsPrefix("foaf",FOAF_NS);
		}
	
	/** reduce the size of the RDF/XML ontology to be loaded */
	private static class OWLFilter implements EventFilter {
		private int in_reject_depth = 0;
		
		private boolean acceptName(final QName qName) {
			final String lcl = qName.getLocalPart();
			if(lcl.equals("Ontology")) return false;
			if(lcl.equals("AnnotationProperty")) return false;
			if(lcl.equals("Axiom")) return false;
			if(lcl.equals("hasDbXref")) return false;
			if(lcl.equals("hasExactSynonym")) return false;
			if(lcl.equals("hasOBONamespace")) return false;
			if(lcl.equals("inSubset")) return false;
			return true;
			}
		@Override
		public boolean accept(final XMLEvent event) {
			if(event.isProcessingInstruction()) return false;
			if(event.isStartElement()) {
				if(!acceptName( event.asStartElement().getName()) ) {
					in_reject_depth++;
					}
				//System.err.println("<"+event.asStartElement().getName().getLocalPart()+"> Returning "+(in_reject_depth==0)+" and after will be "+in_reject_depth);
				return in_reject_depth==0;
				}
			else if(event.isEndElement()) {
				boolean curr  = in_reject_depth==0;
				if( !acceptName(event.asEndElement().getName())) in_reject_depth--;
				//System.err.println("</"+event.asEndElement().getName().getLocalPart()+"> Returning "+curr+" and after will be "+in_reject_depth);
				return curr;
				}
			//System.err.println("("+event.getEventType()+") in_reject_depth "+in_reject_depth);
			return in_reject_depth==0;
			}
		}
	
	
	private abstract static class AbstractEntity {
		private final Resource rsrc;
		protected AbstractEntity(final Resource rsrc) {
			this.rsrc=rsrc;
			}
		public abstract Model getModel();
		
		Resource getResource() {
			return this.rsrc;
			}
		
		@Override
		public boolean equals(final Object obj) {
			if(obj==this) return true;
			if(obj==null || !obj.getClass().equals(this.getClass())) return false;
			return this.rsrc.equals(AbstractEntity.class.cast(obj).rsrc);
			}
		
		@Override
		public final int hashCode() {
			return rsrc.hashCode();
			}
		@Override
		public String toString() {
			return rsrc.toString();
			}
		}
	
	/** class loading  an ONTOLOGY as OWL */
	private class OWLOntology implements AutoCloseable {
		class Term extends AbstractEntity {
			Term(final Resource rsrc) {
				super(rsrc);
				}
			@Override
			public Model getModel() {
				return OWLOntology.this.ontModel;
				}
			
			String getLabel2() {
				final Optional<String> lbl=getLabel();
				final Optional<String> id=getLabel();
				if(id.isPresent() && lbl.isPresent() && id.get().equals(lbl.get())) return id.get();
				return (lbl.isPresent()?lbl.get():"")+"["+id.orElse(getResource().toString())+"]";
				}
			Optional<String> getId() {
				return JenaUtils.stream(getModel().listStatements(getResource(),OWLOntology.this.id_property,(RDFNode)null)).
						map(STMT->STMT.getObject().asLiteral().getString()).
						findFirst();
				}
			Optional<String> getLabel() {
				return JenaUtils.stream(getModel().listStatements(getResource(),RDFS.label,(RDFNode)null)).
						map(STMT->STMT.getObject().asLiteral().getString()).
						findFirst();
				}
			Set<Resource> getDescendantsAsResource() {
				return _descendants(getResource(),new HashSet<>()).
						stream().
						collect(Collectors.toSet());
				}
			
			private Set<Resource> _descendants(final Resource currsrc,final Set<Resource> set) {
				if(set.contains(currsrc)) return set;
				set.add(currsrc);
				JenaUtils.stream(getModel().listResourcesWithProperty(RDFS.subClassOf,currsrc)).
					forEach(R->_descendants(R, set));
				return set;
				}
			
			
			Set<Term> getDescendants() {
				return getDescendantsAsResource().
						stream().
						map(S->new Term(S)).collect(Collectors.toSet());
				}
			boolean isDeprecated() {
				final Literal vrai = getModel().createTypedLiteral(true);
				return 	getModel().contains(getResource(),ResourceFactory.createProperty(OWL.NS+"deprecated"),vrai);
				}
			
			}
		
		
		
		
		
		private final Model ontModel = ModelFactory.createDefaultModel();
		/** how to find term id */
		private final Property id_property;
		private final Path path;
		OWLOntology(final Path path,Property id_property) throws XMLStreamException, IOException, SAXException{
			LOG.info("Loading "+path+"...");
			this.path = path;
			this.id_property = id_property;
			final XMLInputFactory inputFactory = XMLInputFactory.newInstance();
			inputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, true);
			try(Reader r= IOUtils.openFileForReader(path)) {
				final XMLEventReader reader0 = inputFactory.createXMLEventReader(r);
				final XMLEventReader reader = inputFactory.createFilteredReader(reader0, new OWLFilter());
			    StAX2Model.read(reader,this.ontModel,"file://"+path.toString()); 
				reader0.close();
				}
			LOG.info("N="+this.ontModel.size());
			}
		
		Optional<Term> findTerm(final Resource rsrc) {
			return  JenaUtils.stream(this.ontModel.listStatements(rsrc,RDF.type,OWL.Class)).
					map(ST->new Term(ST.getSubject())).
					findAny();
			}
		
		private Set<Term> getClasses() {
			return JenaUtils.stream(this.ontModel.listSubjectsWithProperty(RDF.type,OWL.Class)).
					map(R->new Term(R)).
					collect(Collectors.toSet());
			}
		
		
		void reportError(PrintWriter pw, final Resource rsrc) {
			Optional<Term> t= this.findTerm(rsrc);
			if(!t.isPresent()) {
				pw.println("<"+rsrc+"> was used but not declared in ontology "+this.path);
				}
			else if(t.get().isDeprecated()) {
				pw.println("<"+rsrc+"> was used but it is deprecated in "+this.path);
				}
			}
		@Override
		public void close() {
			this.ontModel.close();
			}
		}
	
	/** base class for sample or group */
	private abstract class AbstractRecordEntity extends AbstractEntity {
		protected AbstractRecordEntity(final Resource rsrc) {
			super(rsrc);
			}
		public final Model getModel() {
			return SamplesRDF.this.model;
			}
		public String getId() {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),PROP_ID)).
					map(O->O.asLiteral().getString()).
					collect(CollectorsUtils.one());
			}
		public String getLabel() {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),RDFS.label)).
					map(O->O.asLiteral().getString()).
					collect(CollectorsUtils.optional()).
					orElse(getId());
			}
		public String getDescription() {
			return StringUtils.ifBlank(JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),RDFS.comment)).
					map(O->O.asLiteral().getString()).
					map(S->StringUtils.normalizeSpaces(S)).
					collect(Collectors.joining(". ")),".");
			}
		}

	/** wrapper for a sample in the model */
	private class Sample extends AbstractRecordEntity {
		protected Sample(final Resource rsrc) {
			super(rsrc);
			}
		
		/** get ID and alias */
		Set<String> getNames() {
			final Set<String> L = new TreeSet<>();
			JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),PROP_alias)).
				filter(O->O.isLiteral()).
				map(O->O.asLiteral().getString()).
				forEach(S->L.add(S));
			L.add(this.getId());
			return L;
			}
		
		Optional<String> getFamily() {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),PROP_family)).
					map(O->O.asLiteral().getString()).
					collect(CollectorsUtils.optional());
			}

		private Optional<Sample> _getParent(final Property link) {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),link)).
					map(O->O.asResource()).
					collect(CollectorsUtils.optional()).
					map(X->new Sample(X));
			}
		OptionalInt getBirthYear() {
			Optional<Integer> x= JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),PROP_birthYear)).
					map(O->O.asLiteral().getInt()).
					collect(CollectorsUtils.optional())
					;
			return x.isPresent()?OptionalInt.of(x.get()):OptionalInt.empty();
			}
		
		Optional<Sample> getFather() {
			return _getParent(PROP_father);
			}
		Optional<Sample> getMother() {
			return _getParent(PROP_mother);
			}
		Optional<Resource> getSex() {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),PROP_sex)).
					filter(O->O.isResource()).
					map(O->O.asResource()).
					collect(CollectorsUtils.optional());
			}
		boolean isMale() {
			final Optional<Resource> sex=getSex();
			return sex.isPresent() && sex.get().equals(IS_MALE);
			}
		boolean isFemale() {
			final Optional<Resource> sex=getSex();
			return sex.isPresent() && sex.get().equals(IS_FEMALE);
			}
		Set<OWLOntology.Term> getTermsInOntology(final OWLOntology ontology,final Property ontology_link) {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),ontology_link)).
					map(O->O.asResource()).
					map(H->ontology.findTerm(H)).
					filter(H->H.isPresent()).
					flatMap(H->H.get().getDescendants().stream()).
					collect(Collectors.toSet())
					;
			}
		Set<Group> getGroups() {
			return JenaUtils.stream(getModel().listSubjectsWithProperty(PROP_contains_sample,getResource())).
					map(R->new Group(R.asResource())).
					collect(Collectors.toSet());
			}
		}

	/** wrapper for a group in the model */
	private class Group extends AbstractRecordEntity {
		protected Group(final Resource rsrc) {
			super(rsrc);
			}
		
		private Set<Group> _getParentGroupAndSelf() {
			return _recursive(this,new HashSet<Group>());
			}
		private Set<Group> _recursive(final Group curr,Set<Group> visited) {
			if(visited.contains(curr)) return visited;
			visited.add(curr);
			JenaUtils.stream(getModel().listObjectsOfProperty(curr.getResource(),PROP_group)).
					map(R->new Group(R.asResource())).
					forEach(X->_recursive(X,visited));
			return visited;
			}
		
		private Set<Sample> _getDeclaredSamples() {
			return JenaUtils.stream(getModel().listObjectsOfProperty(getResource(),PROP_contains_sample)).
					map(R->new Sample(R.asResource())).
					collect(Collectors.toSet());
			}
		Set<Sample> getSamples() {
			return _getParentGroupAndSelf().stream().
					flatMap(G->G._getDeclaredSamples().stream()).
					collect(Collectors.toSet());
			}
		}
	
	/** create entity with id and insert it into model */
	private Resource createEntity(final String id, final Resource type,boolean new_instance) {
		if(id.equals("0") || id.equals(".")) throw new IllegalArgumentException("unsupported id. Cannot be the string '0' or '.'");
		final Resource rsrc = this.model.createResource(NS+"db/" + id.toUpperCase());
		final Optional<Resource> otherType = JenaUtils.stream(model.listStatements(rsrc, RDF.type, rsrc)).
			map(R->R.getObject().asResource()).
			filter(S->!S.equals(TYPE_DECLARED)).
			filter(S->!S.equals(type)).
			findAny();
		if(otherType.isPresent()) {
			throw new IllegalArgumentException(rsrc.toString()+" is already defined as <"+otherType+">. Cannot set type to <"+type+">");
			}
		if(new_instance) {
			 if(model.contains(rsrc, RDF.type, TYPE_DECLARED)) {
				throw new IllegalArgumentException("instance "+rsrc.toString()+" was already defined");
			 	}
			this.model.add(rsrc,RDF.type,TYPE_DECLARED);
			}
		
		this.model.add(rsrc,RDF.type,type);
		this.model.add(rsrc,PROP_ID, id.toUpperCase());
		 return rsrc;
		}
		
	private Resource createHPO(String id0) {
		final String[] tokens = id0.trim().split("[ \t;,]+");
		String id=tokens[0].toUpperCase();
		if(Character.isDigit(id.charAt(0))) id="HP:"+id;
		final Pattern regex = Pattern.compile("HP:[0-9]{7}");
		if(!regex.matcher(id).matches()) throw new IllegalArgumentException("not valid HPO pattern "+regex.pattern()+" for "+id0);
		return ResourceFactory.createResource("http://purl.obolibrary.org/obo/"+id.replace(":", "_"));

		}
	
	private Resource createDOID(String id0) {
		final String[] tokens = id0.trim().split("[ \t;,]+");
		String id=tokens[0].toUpperCase();
		if(Character.isDigit(id.charAt(0))) id="DOID:"+id;
		final Pattern regex = Pattern.compile("DOID:[0-9]{7}");
		if(!regex.matcher(id).matches()) throw new IllegalArgumentException("not valid DOID pattern "+regex.pattern()+" for "+id0);
		return ResourceFactory.createResource("http://purl.obolibrary.org/obo/"+id.replace(":", "_"));
		}
	
	private Resource createSex(String id) {
		id = id.toLowerCase();
		if(id.equals("male")) return IS_MALE;
		if(id.equals("man")) return IS_MALE;
		if(id.equals("female")) return IS_FEMALE;
		if(id.equals("woman")) return IS_FEMALE;
		throw new IllegalArgumentException("not a valid SEX declaration "+id);
		}
	
	

	
	private Map.Entry<String,String> split(final String s) {
		int i= s.indexOf(":");
		if(i==-1) i=s.indexOf("=");
		if(i==-1) throw new IllegalArgumentException("Cannot find key/value delimiter in " + s);
		final String key= s.substring(0,i).trim().toLowerCase();
		if(key.isEmpty()) throw new IllegalArgumentException("Empty key in " + s);
		final String value = s.substring(i+1).trim();
		if(value.isEmpty()) throw new IllegalArgumentException("Empty value in " + s);
		return new AbstractMap.SimpleEntry<>(key,value);
		}
	
	
	private void scan(BufferedReader br,final OWLOntology snomedEthnic) throws IOException {
		final List<String> lines = new ArrayList<>();
		for(;;) {
			String line = br.readLine();
			if(StringUtil.isBlank(line)) {
				if(!lines.isEmpty()) {
					final Resource type = 
							lines.stream().
							map(L->split(L)).
							filter(P->P.getKey().equals("type")).
							map(KV->{
								final String v = KV.getValue().toLowerCase();
								if(v.equals("sample")) return TYPE_SAMPLE;
								if(v.equals("group")) return TYPE_GROUP;
								throw new IllegalArgumentException("unsupported type :"+v+" must be sample or group");
								}).
							collect(CollectorsUtils.one("expect one and only one 'type' per record."))
							;
					
					final String id = 
							lines.stream().
							map(L->split(L)).
							filter(P->P.getKey().equalsIgnoreCase("id")).
							map(KV->KV.getValue().toUpperCase()).
							collect(CollectorsUtils.one("expect one and only one 'id' per record."))
							;
					
					final Resource subject = createEntity(id, type, true);
					
					lines.stream().
						map(L->split(L)).
						map(KV->{
							final String key = KV.getKey();
							if(key.equals("type")) return null;
							if(key.equals("id")) return null;
							if(key.equals("label") || key.equals("name")) return  new AbstractMap.SimpleEntry<Property,RDFNode>(RDFS.label,this.model.createLiteral(KV.getValue()));
							if(key.equals("comment") || key.equals("description")) return  new AbstractMap.SimpleEntry<Property,RDFNode>(RDFS.comment,this.model.createLiteral(KV.getValue()));
							if(key.equals("hpo") || key.equals("hp")) return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_HPO,createHPO(KV.getValue()));
							if(key.equals("doid")) return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_DOID,createDOID(KV.getValue()));

							if(type.equals(TYPE_SAMPLE)) {
								if(key.equals("family")) {
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_family,createEntity(KV.getValue(),TYPE_SAMPLE,false));
									}
								if(key.equals("father")) {
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_father,createEntity(KV.getValue(),TYPE_SAMPLE,false));
									}
								if(key.equals("mother")) {
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_mother,createEntity(KV.getValue(),TYPE_SAMPLE,false));
									}
								if(key.equals("sex")) {
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_sex,createSex(KV.getValue()));
									}
								if(key.equals("birth") || key.equals("birthyear")) {
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_birthYear,model.createTypedLiteral(Integer.parseInt(KV.getValue())));
									}
								if(key.equals("alias")) {
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_alias,model.createLiteral(KV.getValue()));
									}
								if(key.equals("pop") || key.equals("population")) {
									final String pop=StringUtils.normalizeSpaces(KV.getValue());
									if(snomedEthnic!=null) {
										Optional<Resource> optPop = JenaUtils.stream(snomedEthnic.ontModel.listSubjectsWithProperty(RDF.type, OWL.Class)).
											filter(S->JenaUtils.stream(snomedEthnic.ontModel.listObjectsOfProperty(S,RDFS.label)).anyMatch(L->L.asLiteral().getString().equalsIgnoreCase(pop))).
											findFirst();
										
										if(!optPop.isPresent())  optPop = JenaUtils.stream(snomedEthnic.ontModel.listSubjectsWithProperty(RDF.type, OWL.Class)).
												filter(S->JenaUtils.stream(snomedEthnic.ontModel.listObjectsOfProperty(S,SNOMEDID)).anyMatch(L->L.asLiteral().getString().equalsIgnoreCase(pop))).
												findFirst();

										if(optPop.isPresent())  return new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_POPULATION,optPop.get());
										LOG.warn("Cannot find snomed: \""+pop+"\"");
										}
								
									return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_POPULATION,model.createResource(SNOMED_NS + pop.replaceAll("[ ]+","_")));
									}
								}
							if(type.equals(TYPE_GROUP)) {
								if(key.equals("sample")) return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_contains_sample,createEntity(KV.getValue(),TYPE_SAMPLE,false));
								if(key.equals("group") || key.equals("extends")) return  new AbstractMap.SimpleEntry<Property,RDFNode>(PROP_group,createEntity(KV.getValue(),TYPE_GROUP,false));
								}		
							throw new IllegalArgumentException("Cannot handle "+key+" for type <"+type+">.");
							}).
						filter(KV->KV!=null).
						forEach(KV->model.add(subject,KV.getKey(),KV.getValue()));
						;
					
					}
				if(line==null) break;
				lines.clear();
				continue;
				}
			if(line.startsWith("#")) continue;
			lines.add(line);
			}
		}
	
	private void saveOntology(final OWLOntology ontology,final String fileSuffix,final Property link_property) throws IOException{
		if(ontology==null) return;
		final Path path = Paths.get(this.outputBase + fileSuffix);
		LOG.info("Saving "+path+"...");
		try(PrintWriter pw = IOUtils.openPathForPrintWriter(path)) {
			pw.println("#uri id label samples".replace(' ', '\t'));
			for(OWLOntology.Term term: ontology.getClasses()) {
				final Set<Sample> samples = JenaUtils.stream(this.model.listSubjectsWithProperty(RDF.type,TYPE_SAMPLE)).
						map(SN->new Sample(SN)).
						filter(R->R.getTermsInOntology(ontology, link_property).contains(term)).
						collect(Collectors.toSet());
				if(samples.isEmpty()) continue;
				pw.print(term.getResource().toString());
				pw.print("\t");
				pw.print(term.getId().orElse("."));
				pw.print("\t");
				pw.print(term.getLabel().orElse("."));
				pw.print("\t");
				pw.print(samples.stream().map(S->S.getId()).sorted().collect(Collectors.joining(",")));
				pw.println();
				}
			pw.flush();
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			final OWLOntology hpoOntology = (this.hpoFilename==null?null:new OWLOntology(this.hpoFilename,OBOINOWBL_id));
			final OWLOntology doidOntology = (this.doidFilename==null?null:new OWLOntology(this.doidFilename,OBOINOWBL_id));
			final OWLOntology snomedEthnic = (this.snomedEthnicFilename==null?null:new OWLOntology(this.snomedEthnicFilename,SNOMEDID));
			
			
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				try(BufferedReader br = IOUtils.openStreamForBufferedReader(stdin())) {
					scan(br,snomedEthnic);
					}
				}
			else
				{
				for(Path path: paths) {
					try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
						scan(br,snomedEthnic);
						}
					}
				}
			final List<Group> all_groups = JenaUtils.stream(model.listSubjectsWithProperty(RDF.type,TYPE_GROUP)).
					map(R->new Group(R)).
					collect(Collectors.toSet()).
					stream().
					sorted((A,B)->A.getId().compareTo(B.getId())).
					collect(Collectors.toList());
			
			final List<Sample> all_samples = JenaUtils.stream(model.listSubjectsWithProperty(RDF.type,TYPE_SAMPLE)).
					map(R->new Sample(R)).
					collect(Collectors.toSet()).
					stream().
					sorted((A,B)->A.getId().compareTo(B.getId())).
					collect(Collectors.toList());
			
			//finalize transfer DOID and HPO from group to samples
			for(Group G:all_groups) {
				final Set<Resource> hpos = JenaUtils.stream(model.listObjectsOfProperty(G.getResource(),PROP_HPO)).
						filter(O->O.isResource()).
						map(O->O.asResource()).
						collect(Collectors.toSet());
				final Set<Resource> doids = JenaUtils.stream(model.listObjectsOfProperty(G.getResource(),PROP_DOID)).
						filter(O->O.isResource()).
						map(O->O.asResource()).
						collect(Collectors.toSet());
				for(Sample sn:G.getSamples()) {
					for(Resource hpo:hpos) {
						this.model.add(sn.getResource(),PROP_HPO,hpo);
						this.model.remove(G.getResource(),PROP_HPO,hpo);
						}
					for(Resource doid:doids) {
						this.model.add(sn.getResource(),PROP_DOID,doid);
						this.model.remove(G.getResource(),PROP_DOID,doid);
						}
					}
				}
			
			
			
			/** save errors */
			try(PrintWriter pw = IOUtils.openPathForPrintWriter(Paths.get(this.outputBase + "_errors.txt"))) {
				//post validate
				JenaUtils.stream(model.listStatements(null,RDF.type,(RDFNode)null)).
					filter(S->S.getObject().equals(TYPE_SAMPLE) || S.getObject().equals(TYPE_GROUP)).
					filter(STMT->!model.contains(STMT.getSubject(),RDF.type,TYPE_DECLARED)).
					forEach(STMT->{
						pw.println("<"+STMT.getSubject()+"> rdf:type <"+STMT.getObject()+"> was used but not declared");
					});
				
				if(hpoOntology!=null) {
					JenaUtils.stream(model.listObjectsOfProperty(PROP_HPO)).
						map(O->O.asResource()).
						forEach(R->hpoOntology.reportError(pw,R));
					}
				if(doidOntology!=null) {
						JenaUtils.stream(model.listObjectsOfProperty(PROP_DOID)).
						map(O->O.asResource()).
						forEach(R->doidOntology.reportError(pw,R));
					}
				if(snomedEthnic!=null) {
					JenaUtils.stream(model.listObjectsOfProperty(PROP_POPULATION)).
					map(O->O.asResource()).
					forEach(R->snomedEthnic.reportError(pw,R));
					}
				for(int x=0;x+1 < all_samples.size();x++) {
					final Set<String> set1 = all_samples.get(x).getNames();
					for(int y=x+1; y < all_samples.size();y++) {
						final Set<String> set2 = all_samples.get(y).getNames();
						set2.retainAll(set1);
						if(!set2.isEmpty()) {
							pw.println("ERROR sample "+all_samples.get(x)+" and "+all_samples.get(y)+" share the same name/alias:"+String.join(",",set2));
						}
					}

				}
				
				for(Group G:all_groups) {
					if(G.getSamples().isEmpty()) {
						pw.println("group "+G.getId()+" is empty.");
						}
					}
				for(Sample R:all_samples) {
					Optional<String> fam = R.getFamily();
					OptionalInt year = R.getBirthYear();
					for(int side=0;side<2;side++) {
						final Sample parent = (side==0?R.getFather():R.getMother()).orElse(null);
						if(parent!=null) {
							if(!(side==0?parent.isMale():parent.isFemale())) {
								pw.println(parent.getId()+" declared as "+(side==0?"father":"mother")+" of "+R.getId()+". But sex is not "+(side==0?"male":"female"));
								}
							Optional<String> famp = parent.getFamily();
							OptionalInt fyear = parent.getBirthYear();
							if(fam.isPresent() && famp.isPresent() && !fam.get().equals(famp.get())) {
								pw.println(parent.getId()+" don't have the same family than"+R+" ("+famp+" vs "+ fam+")");
								}
							else if(!fam.isPresent() && !famp.isPresent() ) {
								// nothing
								}
							else
								{
								pw.println(parent.getId()+" and "+R+" don't have both family specified ("+famp+" vs "+ fam+")");
								}
							if(year.isPresent() && fyear.isPresent() && year.getAsInt()>=fyear.getAsInt()) {
								pw.println(parent.getId()+" is younder than is child "+R+" ("+year+" vs "+ fyear+")");
								}
							}
						}
					} // samples
				} //errors
				
			/* at this point cleanup the data */
			this.model.removeAll(null, RDF.type, TYPE_DECLARED);
			try(PrintWriter pw = IOUtils.openPathForPrintWriter(Paths.get(this.outputBase + "_model.rdf"))) {
				model.write(pw,"RDF/XML");
				pw.flush();
				}
			
			/** save DOID/HPO output file */
			saveOntology(doidOntology, "_doid.tsv", PROP_DOID);
			saveOntology(hpoOntology, "_hpo.tsv", PROP_HPO);
			saveOntology(snomedEthnic, "_populations.tsv", PROP_POPULATION);
			
			if(!all_groups.isEmpty()) {
				LOG.info("saving groups");
				try(PrintWriter pw = IOUtils.openPathForPrintWriter(Paths.get(this.outputBase + "_groups.tsv"))) {
					pw.println("#id label desc samples.count male.count female.count other.count groups".replace(' ', '\t'));
					for(Group R: all_groups) {
						pw.print(R.getId());
						pw.print("\t");
						pw.print(R.getLabel());
						pw.print("\t");
						pw.print(R.getDescription());
						pw.print("\t");
						pw.print(R.getSamples().stream().map(G->G.getId()).count());
						pw.print("\t");
						pw.print(R.getSamples().stream().filter(G->G.isMale()).count());
						pw.print("\t");
						pw.print(R.getSamples().stream().filter(G->G.isFemale()).count());
						pw.print("\t");
						pw.print(R.getSamples().stream().filter(G->!(G.isFemale() || G.isMale())).count());
						pw.print("\t");
						pw.print(StringUtils.ifBlank(R.getSamples().stream().map(G->G.getId()).sorted().collect(Collectors.joining(";")),"."));
						pw.println();
						};
					pw.flush();
					}
				} /* end group */
			
			if(!all_samples.isEmpty()) {
				LOG.info("saving samples");
				try(PrintWriter pw = IOUtils.openPathForPrintWriter(Paths.get(this.outputBase + "_samples.tsv"))) {
					pw.println("#id label desc family aliases birthYear father mother sex HPO DOID POP groups".replace(' ', '\t'));
					for(Sample R: all_samples) {
							pw.print(R.getId());
							pw.print("\t");
							pw.print(R.getLabel());
							pw.print("\t");
							pw.print(R.getDescription());
							pw.print("\t");
							pw.print(R.getFamily().orElse("."));
							pw.print("\t");
							pw.print(StringUtils.ifBlank(R.getNames().stream().filter(S->!S.equals(R.getId())).sorted().collect(Collectors.joining(",")),"."));
							pw.print("\t");
							if(R.getBirthYear().isPresent()) {
								pw.print(R.getBirthYear().getAsInt());
								}
							else
								{
								pw.print(".");
								}
							pw.print("\t");
							Optional<Sample> parent = R.getFather();
							pw.print(parent.isPresent()?parent.get().getId():".");
							pw.print("\t");
							parent = R.getMother();
							pw.print(parent.isPresent()?parent.get().getId():".");
							pw.print("\t");
							if(R.isMale()) {
								pw.print("male");
								}
							else if(R.isFemale()) {
								pw.print("female");
								}
							else
								{
								pw.print(".");
								}
							pw.print("\t");
							if(hpoOntology!=null) {
								pw.print(StringUtils.ifBlank(R.getTermsInOntology(hpoOntology, PROP_HPO).stream().
									map(T->T.getLabel2()).
									collect(Collectors.joining("; "))
									,"."));
								}
							else
								{
								pw.print(".");
								}
							pw.print("\t");
							if(doidOntology!=null) {
								pw.print(StringUtils.ifBlank(R.getTermsInOntology(doidOntology, PROP_DOID).stream().
										map(T->T.getLabel2()).
										collect(Collectors.joining("; "))
										,"."));
								}
							else
								{
								pw.print(".");
								}
							pw.print("\t");
							if(snomedEthnic!=null){
								pw.print(StringUtils.ifBlank(R.getTermsInOntology(snomedEthnic, PROP_POPULATION).stream().
										map(T->T.getLabel2()).
										collect(Collectors.joining("; "))
										,"."));
								}
							else
								{
								pw.print(".");
								}
							pw.print("\t");
							pw.print(StringUtils.ifBlank(R.getGroups().stream().map(G->G.getId()).sorted().collect(Collectors.joining(";")),"."));
							pw.println();
						}
					pw.flush();
					}
				}//end samples
			
				if(!all_samples.isEmpty()) {
				/** save pedigree */
				LOG.info("saving pedigree");
					try(PrintWriter pw = IOUtils.openPathForPrintWriter(Paths.get(this.outputBase + "_pedigree.txt"))) {
						for(Sample R: all_samples) {
							final String id = R.getId();
							pw.print(R.getFamily().orElse(id));
							pw.print("\t");
							pw.print(id);
							pw.print("\t");
							Optional<Sample> father = R.getFather();
							pw.print(father.isPresent()?father.get().getId():"0");
							pw.print("\t");
							Optional<Sample> mother = R.getMother();
							pw.print(mother.isPresent()?mother.get().getId():"0");
							pw.print("\t");
							if(R.isMale()) {
								pw.print("1");
								}
							else if(R.isFemale()) {
								pw.print("2");
								}
							else
								{
								pw.print("0");
								}	
							pw.println();
								}
						pw.flush();
						}
					} /* save pedigree */
				
				if(!all_samples.isEmpty() && !all_groups.isEmpty()) {
					/** save grop 2 sample */
					LOG.info("saving group to samples");
					try(PrintWriter pw = IOUtils.openPathForPrintWriter(Paths.get(this.outputBase + "_group2sample.txt"))) {
						pw.println("#group sample".replace(' ', '\t'));
						JenaUtils.stream(model.listSubjectsWithProperty(RDF.type,TYPE_GROUP)).
							map(R->new Group(R)).
							forEach(G->{
								for(Sample sn:G.getSamples()) {
									pw.print(G.getId());
									pw.print("\t");
									pw.print(sn.getId());
									pw.println();
									}
								});
						pw.flush();
						}
					}/* end save group/sample */
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			this.model.close();
			}
		}
	
	public static void main(String[] args) {
		new SamplesRDF().instanceMainWithExit(args);

	}

}
