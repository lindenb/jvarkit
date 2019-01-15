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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;


public class Pedigree
	{
	private static final Logger LOG = Logger.build(Pedigree.class).make();
	public static final String OPT_DESCRIPTION="A pedigree is a text file delimited with tabs. No header. Columns are (1) Family (2) Individual-ID (3) Father Id or '0' (4) Mother Id or '0' (5) Sex : 1 male/2 female / 0 unknown (6) Status : 0 unaffected, 1 affected,-9 unknown ";
	private final Map<String,FamilyImpl > families=new TreeMap<String, Pedigree.FamilyImpl>();
	
	public enum Status
		{
		missing,unaffected,affected;
		public int intValue()
			{
			switch(this)
				{
				case missing: return -9;
				case unaffected: return 0;
				case affected: return 1;
				default:throw new IllegalStateException();
				}
			}
		}
	
	public enum Sex
		{
		male,female,unknown;
		public int intValue()
			{
			switch(this)
				{
				case male: return 1;
				case female: return 2;
				case unknown: return 0;
				default:throw new IllegalStateException();
				}
			}
		}
	
	public interface Family extends Comparable<Family>
		{
		public String getId();
		public Person getPersonById(String s);
		public java.util.Collection<? extends Person> getIndividuals();
		public Family validate() throws IllegalStateException;
		@Override
		default int compareTo(final Family o) {
			return this.getId().compareTo(o.getId());
			}
		}
	
	private class FamilyImpl
		implements Family
		{
		private String id;
		private Map<String,PersonImpl> individuals=new TreeMap<String,PersonImpl>();
		
		@Override
		public String getId()
			{
			return this.id;
			}
		@Override
		public Person getPersonById(String s)
			{
			return individuals.get(s);
			}
		@Override
		public java.util.Collection<? extends Person> getIndividuals()
			{
			return this.individuals.values();
			}
		@Override
		public int hashCode() {
			return id.hashCode();
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null || getClass() != obj.getClass()) {
				return false;
			}
			
			final FamilyImpl other = (FamilyImpl) obj;
			return id.equals(other.id);
			}
		@Override
		public String toString() {
			return this.id;
			}
		
		@Override
		public Family validate() throws IllegalStateException {
			for(final PersonImpl p: this.individuals.values()) {
				p.validate();
				}
			return this;
			}
		
		}
	
	public interface Person
		extends Comparable<Person>
		{
		public String getId();
		/** the ID of this individual is unique in the pedigree. It can be associated without ambiguity to a sample in a VCF header */
		public boolean hasUniqId() ;
		public Family getFamily();
		public boolean hasFather();
		public Person getFather();
		public boolean hasMother();
		public Person getMother();
		/** returns true if father or mother is present in pedigree */
		public default boolean hasAtLeastOneParent()
			{
			return hasFather() || hasMother();
			}
		
		/** returns the father id or null if there is no father */
		public default String getFatherId() { 
			final Person p = getFather();
			return p==null?null:p.getId();
			}
		public default String getMotherId() { 
			final Person p = getMother();
			return p==null?null:p.getId();
			}
		
		public Sex getSex();
		public boolean isMale();
		public boolean isFemale();
		public Status getStatus();
		/** return i-th parent 0=father 1=mother */
		public Person getParent( int zeroOrOne);
		public Person validate() throws IllegalStateException;
		
		public boolean isAffected();
		public boolean isUnaffected();
		
		@Override
		default int compareTo(final Person o) {
			int i= getFamily().compareTo(o.getFamily());
			if(i!=0) return i;
			return this.getId().compareTo(o.getId());
			}
		
		}
	
	
	private class PersonImpl
	implements Person
		{
		FamilyImpl family;
		String id;
		String fatherId=null;
		String motherId=null;
		Sex sex=Sex.unknown;
		Status status=Status.unaffected;
		
		@Override
		public Family getFamily()
			{
			return this.family;
			}
		
		@Override
		public String getId()
			{
			return this.id;
			}
		
		@Override
		public boolean isMale() { return Sex.male.equals(this.getSex());}
		@Override
		public boolean isFemale() { return Sex.female.equals(this.getSex());}

		@Override
		public boolean isAffected() { return Status.affected.equals(this.getStatus());}
		@Override
		public boolean isUnaffected() { return Status.unaffected.equals(this.getStatus());}

		
		
		private Person getParent(final String s)
			{
			if(s==null || s.isEmpty() || s.equals("0")) return null;
			return getFamily().getPersonById(s);
			}
		
		public Person getParent( int zeroOrOne) {
			switch(zeroOrOne) {
			case 0: return getFather();
			case 1: return getMother();
			default: throw new IllegalArgumentException("0 or 1 but got "+zeroOrOne);
			}
		}
		
		private boolean hasParent(final String s)
			{
			return !(s==null || s.isEmpty() || s.equals("0"));
			}
		
		@Override
		public boolean hasFather() {
			return hasParent(this.fatherId);
		}
		
		@Override
		public boolean hasMother() {
			return hasParent(this.motherId);
		}
		
		public Person getFather()
			{
			return getParent(fatherId);
			}
		
		public Person getMother()
			{
			return getParent(motherId);
			}
		
		@Override
		public Status getStatus()
			{
			return status;
			}
		
		@Override
		public Sex getSex()
			{
			return sex;
			}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + family.hashCode();
			result = prime * result + id.hashCode();
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj) return true;
			if (obj == null || getClass() != obj.getClass()) {
				return false;
			}
			final PersonImpl p =(PersonImpl)obj;
			return	this.family.equals(p.family) &&
					this.id.equals(p.id);
			}
		
		@Override
		public boolean hasUniqId() 
			{
			return (Pedigree.this.families.values().
					stream().
					flatMap(F->F.getIndividuals().stream())).
					filter(P->getId().equals(P.getId())).count() ==1L;
			}

		@Override
		public Person validate() throws IllegalStateException {
			if(this.fatherId!=null && !this.fatherId.equals("0")) {
				final PersonImpl parent = this.family.individuals.get(this.fatherId);
				if(parent==null) throw new IllegalStateException(
					"Individual "+this.toString()+" has father "+this.fatherId+" "+
					"but he is missing in family."
					);
				if(parent.sex == Sex.female) throw new IllegalStateException(
					"Individual "+this.toString()+" has father "+this.fatherId+" "+
					"but he is declared as a woman."
					);
			}
			if(this.motherId!=null && !this.motherId.equals("0")) {
				final PersonImpl parent = this.family.individuals.get(this.motherId);
				if(parent==null) throw new IllegalStateException(
					"Individual "+this.toString()+" has mother "+this.motherId+" "+
					"but she is missing in family."
					);
				if(parent.sex == Sex.male) throw new IllegalStateException(
					"Individual "+this.toString()+" has mother "+this.motherId+" "+
					"but she is declared as a man."
					);
			}
			return this;
		}
		
		@Override
		public String toString() {
			return family+":"+this.id;
			}
		}
	
	
	
	private Pedigree()
		{
		
		}
	
	public static Pedigree createEmptyPedigree() {
		return new Pedigree();
	}
	
	public boolean isEmpty() {
		return this.families.isEmpty();
	}
	
	/** returns true if this pedigree has at least one trio (  children has at leat one parent in the pedigree)*/
	public boolean hasTrios()
		{
		return this.getPersons().stream().
				anyMatch(P->P.hasFather() || P.hasMother());
				
		}
	/** return a list of Children in trio is this pedigree . children has at leat one parent in the pedigree*/
	public List<Person> getChildrenInTrios()
		{
		return this.getPersons().stream().
				filter(P->P.hasFather() || P.hasMother()).
				collect(Collectors.toList());
		}
	
	/** utility function for vcf, returns true if all person's ID in the pedigree
	 * are unique (no same ID shared by two families 
	 */
	public boolean verifyPersonsHaveUniqueNames() {
		final Set<String> m = new HashSet<String>();
		for(final Family f:families.values())
			{
			for(final Person p:f.getIndividuals()) {
				if(m.contains(p.getId())) {
					return false;
					}
				m.add(p.getId());
				}
			}
		return true;
		}
	
	/** utility function for vcf, return a Map<person.id,Person>
	 * will throw an illegalState if two individual have the same ID (but ! families )
	 * @return Map<person.id,Person>
	 */
	public Map<String, Person> getPersonsMap() {
		final Map<String, Person> m = new TreeMap<String, Person>();
		for(final Family f:families.values())
			{
			for(final Person p:f.getIndividuals()) {
				final Person prev = m.get(p.getId());
				if(prev!=null) {
					throw new IllegalStateException(
						"Cannot create a Map<String, Person> because "+prev+" and "+p + " share the same ID : "+p.getId()	
						);
					}
				m.put(p.getId(), p);
				}
			}
		return m;
	}
	
	/** validate pedigree */
	public Pedigree validate() throws IllegalStateException {
		for(final Family f: getFamilies()) f.validate();
		return this;
		}

	
	/** get all the families in this pedigree */
	public java.util.Collection<? extends Family> getFamilies()
		{
		return this.families.values();
		}
	
	public Family getFamilyById(final String famId) {
		return this.families.get(famId);
	}
	
	/** get an individual by id, assume individual-ID are unique */
	public Person getPersonById(final String id) {
		for(final Family fam:this.families.values())
			{
			final Person p = fam.getPersonById(id);
			if(p!=null) return p;
			}
		return null;
	}

	/** get an individual by id, assume individual-ID are unique, return null on not found */
	public Person getUniqPersonById(final String id) {
		Person p = null;
		for(final Family fam:this.families.values())
			{
			final Person p2 = fam.getPersonById(id);
			if(p2==null) continue;
			if(p!=null) throw new IllegalStateException("found duplicate sample by id "+p+" and "+p2); 
			p = p2;
			}
		return p;
	}

	
	
	/** get all the individuals in this pedigree */
	public java.util.Set<Person> getPersons()
		{
		final java.util.Set<Person> set = new HashSet<>();
		for(final Family f:families.values())
			{
			set.addAll(f.getIndividuals());
			}
		return set;
		}
	
	/** contains affected individuals */
	public boolean hasAffected()
		{
		return getPersons().stream().
				anyMatch(P->P.getStatus()==Status.affected)
				;
		}
	
	/** get affected individuals */
	public java.util.Set<Person> getAffected()
		{
		return getPersons().stream().filter(P->P.getStatus()==Status.affected).collect(Collectors.toSet());
		}

	/** contains unaffected individuals */
	public boolean hasUnaffected()
		{
		return getPersons().stream().
				anyMatch(P->P.getStatus()==Status.unaffected)
				;
		}
	
	/** get unaffected individuals */
	public java.util.Set<Person> getUnaffected()
		{
		return getPersons().stream().filter(P->P.getStatus()==Status.unaffected).collect(Collectors.toSet());
		}
	
	@Override
	public String toString() {
		return "<<Pedigree>>";
		}
	
	@Deprecated //use Pedigree.Parser.parse
	public static Pedigree readPedigree(final File f) throws IOException
		{
		return newParser().parse(f);
		}	
	@Deprecated //use Pedigree.Parser.parse
	public static Pedigree readPedigree(final BufferedReader r) throws IOException
		{
		return newParser().parse(r);
		}	
	
	public static final String VcfHeaderKey="Sample";
	
	@Deprecated //use Pedigree.Parser.parse
	public static Pedigree readPedigree(final VCFHeader header) {
		return newParser().parse(header);
		}
	@Deprecated //use Pedigree.Parser.parse
	public static Pedigree readPedigree(final Collection<VCFHeaderLine> metadata) {
		return newParser().parse(metadata);
	}
	
	public Set<VCFHeaderLine> toVCFHeaderLines()  {
		final Set<VCFHeaderLine> set = new LinkedHashSet<>();
		for(final Family f:families.values())
			{
			for(final Person p:f.getIndividuals()) {
				final StringBuilder sb=new StringBuilder();
				sb.append("<Family=");
				sb.append(f.getId());
				sb.append(",ID=");
				sb.append(p.getId());
				sb.append(",Father=");
				sb.append(p.getFather()==null?"0":p.getFather().getId());
				sb.append(",Mother=");
				sb.append(p.getMother()==null?"0":p.getMother().getId());
				sb.append(",Sex=");
				sb.append(p.getSex().intValue());
				sb.append(",Status=");
				sb.append(p.getStatus().intValue());
				sb.append(">");
				set.add(new VCFHeaderLine(VcfHeaderKey, sb.toString()));
				}
			}
		
		return set;
		}
	

	
	/** creates a new PEdigree parser */
	public static Parser newParser()
		{
		return new Parser();
		}
	
	/** how to interpret the 'affected' column */
	public enum StatusModel implements Function<String, Status>{
		un0af1() {
			@Override
			public Status apply(String status) {
				if(status==null) return null;
				if(status.equals("1")) return Status.affected;
				else if(status.equals("0")) return Status.unaffected;	
				return null;
				}
		};
		
	};
	
	public static final StatusModel DefaultStatusModel = StatusModel.un0af1;
	
	public static class Parser
		{
		private final Pattern tab = Pattern.compile("[\t]");
		private StatusModel statusModel = Pedigree.DefaultStatusModel;
		private boolean statusRequired=false;
		private boolean acceptErrorsInVcfHeader=false;
		
		public Pedigree parse(final File f) throws IOException
			{
			try(BufferedReader r= IOUtils.openFileForBufferedReading(f)) {
				return this.parse(r);
				}
			}
		public Pedigree parse(final BufferedReader r)throws IOException
			{
			final Pedigree ped = new Pedigree();
			r.lines().forEach(L->{
				if(L.isEmpty() || L.startsWith("#")) return;
				read(ped,tab.split(L));
				});
			return ped;
			}
		
		public Parser statusModel(StatusModel statusModel)
			{
			this.statusModel = statusModel;
			return this;
			}
		public Parser statusIsRequired(boolean statusRequired)
			{
			this.statusRequired = statusRequired;
			return this;
			}
		public Parser acceptErrorsInVcfHeader(boolean acceptErrorsInVcfHeader)
			{
			this.acceptErrorsInVcfHeader = acceptErrorsInVcfHeader;
			return this;
			}
		private void read(final Pedigree ped,final String tokens[])
			{
			final String fam= tokens[0];
			final String indi = tokens[1];
			final String father = tokens[2];
			final String mother = tokens[3];
			final String sex = (tokens.length>4?tokens[4]:"");
			final String status = (tokens.length>5?tokens[5]:"");
			build(ped,fam,indi,father,mother,sex,status);
			}
		
		
		public Pedigree parse(final VCFHeader h)
			{
			return this.parse(h.getMetaDataInInputOrder());
			}
		
		private void build(final Pedigree ped,final String famId,final String indiId,final String fatherId,final String motherId,final String sexxx,final String status)
			{
			FamilyImpl fam=ped.families.get(famId);
			if(fam==null)
				{
				fam=ped.new FamilyImpl();
				fam.id = famId;
				ped.families.put(famId, fam);
				}
			if(fam.getPersonById(indiId)!=null) throw new IllegalArgumentException("duplicate individual: "+String.join(" ; ", famId,indiId,fatherId,motherId,sexxx,status));
			final PersonImpl p=ped.new PersonImpl();
			p.family=fam;
			p.id=indiId;
			p.fatherId=fatherId;
			p.motherId=motherId;
			
			if(sexxx!=null)
				{
				if(sexxx.equals("1")) p.sex=Sex.male;
				else if(sexxx.equals("2")) p.sex=Sex.female;
				}
			
			if(status!=null)
				{
				final Status st= this.statusModel.apply(status);
				if(st!=null ) p.status=st;
				}
			else if(this.statusRequired) {
				throw new IllegalArgumentException("status must be declared");
				}
			
			
			fam.individuals.put(p.id, p);		
			}

		
		/** should be readPedigree(header.getMetaDataInInputOrder()) */
		public Pedigree parse(final Collection<VCFHeaderLine> metadata)
			{
			final java.util.function.Consumer<String> handleError= (S)->{
				if(acceptErrorsInVcfHeader)
					{
					LOG.warn(S);
					}
				else
					{
					throw new IllegalArgumentException(S);
					}
				};
			
			final Pattern comma = Pattern.compile("[,]");
			final Pedigree ped=new Pedigree();
			for(final VCFHeaderLine h:metadata) {
				final String key = h.getKey();
				if(!VcfHeaderKey.equals(key)) continue;
				final String value =h.getValue();
				if(!value.startsWith("<")) {
					LOG.warn("in "+VcfHeaderKey+" value doesn't start with '<' "+value);
					continue;
				}
				if(!value.endsWith(">")) {
					LOG.warn("in "+VcfHeaderKey+" value doesn't end with '>' "+value);
					continue;
				}
				
				String familyId=null;
				String indiId=null;
				String fatherId=null;
				String motherId=null;
				String sexx=null;
				String status=null;

				for(final String t:comma.split(value.substring(1, value.length()-1))) {
					final int eq = t.indexOf("=");
					if(eq==-1) 
						{
						LOG.warn("'=' missing in "+t+" of "+value);
						continue;
						}
					final String left = t.substring(0,eq);
					if(left.equals("Family")) {
						if(familyId!=null) {
							handleError.accept("Family defined twice in " +value);
							continue;
							}
						familyId= t.substring(eq+1).trim();
						}
					else if(left.equals("ID")) {
						if(indiId!=null) {
							handleError.accept("ID defined twice in " +value);
							continue;
							}
						indiId= t.substring(eq+1).trim();
						}
					else if(left.equals("Father")) {
						if(fatherId!=null) {
							handleError.accept("fatherId defined twice in " +value);
							continue;
						}
						fatherId= t.substring(eq+1).trim();
						}
					else if(left.equals("Mother")) {
						if(motherId!=null) {
							handleError.accept("mother defined twice in " +value);
							continue;
						}
						motherId= t.substring(eq+1).trim();
						}
					else if(left.equals("Sex")) {
						if(sexx!=null) {
							handleError.accept("sex defined twice in " +value);
							continue;
							}
						sexx= t.substring(eq+1).trim();
						}
					else if(left.equals("Status")) {
						if(status!=null) {
							handleError.accept("status defined twice in " +value);
							continue;
						}
						status= t.substring(eq+1).trim();
						}
					}
				if(familyId==null) {
					handleError.accept("Family undefined  in " +value);
					continue;
				}
				if(indiId==null) {
					handleError.accept("ID undefined in " +value);
					continue;
				}
				build(ped,familyId,indiId,fatherId,motherId,sexx,status);
				}
			return ped;
			}
		}
	
	/** extract case controls in VCF header injected with VcfInjectPedigree */
	public static class CaseControlExtractor
		{
		private boolean verbose= true;
		
		public void setVerbose(boolean verbose) {
			this.verbose = verbose;
			}
		public boolean isVerbose() {
			return verbose;
		}
		
		/** extract affected/non-affected sample in pedigree that exists in VCF header samples */
		public java.util.Set<com.github.lindenb.jvarkit.util.Pedigree.Person> extract(
				final htsjdk.variant.vcf.VCFHeader header,
				final com.github.lindenb.jvarkit.util.Pedigree pedigree
				)
			{
			if(!pedigree.verifyPersonsHaveUniqueNames()) {
				throw new JvarkitException.PedigreeError("I can't use this pedigree in a VCF because two samples have the same ID.");
			}

			final java.util.Set<String> samplesNames= new java.util.HashSet<>(header.getSampleNamesInOrder());
			final java.util.Set<com.github.lindenb.jvarkit.util.Pedigree.Person> individuals = new java.util.HashSet<>(pedigree.getPersons());
			
			
			final java.util.Iterator<com.github.lindenb.jvarkit.util.Pedigree.Person> iter= individuals.iterator();
			while(iter.hasNext())
			{
				final com.github.lindenb.jvarkit.util.Pedigree.Person person = iter.next();
				if(!(samplesNames.contains(person.getId()) && (person.isAffected() || person.isUnaffected()))) {
					if(isVerbose()) LOG.warn("Ignoring "+person+" because it is not present in VCF header or status is unknown");
					iter.remove();
				}
			}
			
			if(isVerbose())  {
				LOG.info("Individuals :"+individuals.size() +
					" affected :"+individuals.stream().filter(P->P.isAffected()).count() +
					" unaffected :"+individuals.stream().filter(P->P.isUnaffected()).count()
					);
				}
			return java.util.Collections.unmodifiableSet( individuals );
		}
		
		public java.util.Set<com.github.lindenb.jvarkit.util.Pedigree.Person> extract(final htsjdk.variant.vcf.VCFHeader header) {
			final com.github.lindenb.jvarkit.util.Pedigree pedigree = Pedigree.newParser().parse(header);
			if(pedigree.isEmpty())
				{
				throw new JvarkitException.PedigreeError("No pedigree found in header. use VcfInjectPedigree to add it");
				}
			return extract(header,pedigree);			
			}	
		}
	

	}
