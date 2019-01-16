package com.github.lindenb.jvarkit.pedigree;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Pedigree;

import sun.security.util.Length;


public class PedigreeParser {
	
private Function<String, Sex> sexParser = S->{
	if(S.length()==1) {
		switch(S.charAt(0)) {
			case '0': return Sex.unknown;
			case '1': return Sex.male;
			case '2': return Sex.female;
			}
		}
	throw new IllegalArgumentException("bad sex in \""+S+"\", should be 0/1/2.");
	};
	
/** creates an empty Pedigree */
public static Pedigree empty() {
	return new PedigreeImpl();
	}

private static class PedigreeImpl implements Pedigree {
	private final Map<String,Family> id2families = new TreeMap<>();
	@Override
	public Collection<Family> getFamilies() {
		return this.id2families.values();
		}
	public Family getFamilyById(final String id) {
		return this.id2families.get(id);
		}
	}

public Pedigree parse(final BufferedReader br)throws IOException
	{
	final CharSplitter tab=CharSplitter.TAB;
	final PedigreeImpl ped = new PedigreeImpl();
	String line;
	while((line=br.readLine())!=null) {
		if(StringUtils.isBlank(line)) continue;
		if(line.startsWith("#")) continue;
		final String tokens[]= tab.split(line);
		if(tokens.length<4) throw new IllegalArgumentException();
		final String famId= tokens[0];
		final String indiId = tokens[1];
		final String fatherId = tokens[2];
		final String motherId = tokens[3];
		final String sex = (tokens.length>4?tokens[4]:"");
		final String status = (tokens.length>5?tokens[5]:"");
		build(ped, famId, indiId, fatherId, motherId, sex, status);
		}
	return ped;
	}

private void build(final PedigreeImpl ped,final String famId,final String indiId,final String fatherId,final String motherId,final String sexxx,final String status)
	{
	FamilyImpl fam= FamilyImpl.class.cast(ped.id2families.get(famId));
	if(fam==null)
		{
		fam= new FamilyImpl(ped,famId);
		ped.id2families.put(famId, fam);
		}
	if(fam.id2individuals.containsKey(indiId)) throw new IllegalArgumentException("duplicate individual: "+String.join(" ; ", famId,indiId,fatherId,motherId,sexxx,status));
	final PersonImpl p= new PersonImpl(fam,indiId);
	p.fatherId=fatherId;
	p.motherId=motherId;
	
	if(sexxx!=null)
		{
		p.sex = this.sexParser.apply(sexxx);
		}
	
	if(status!=null)
		{
		final Status st= this.statusModel.apply(status);
		if(st!=null ) p.status=st;
		}
	else if(this.statusRequired) {
		throw new IllegalArgumentException("status must be declared");
		}
	
	
	fam.id2individuals.put(p.id, p);		
	}

private static class FamilyImpl implements Family
	{
	final PedigreeImpl ped;
	private final String id;
	private Map<String,SampleImpl> id2individuals=new TreeMap<>();
	private FamilyImpl(final PedigreeImpl ped,final String id) {
		this.ped = ped;
		this.id = id;
	}
	@Override
	public String getId()
		{
		return this.id;
		}
	@Override
	public Sample getSampleById(final String s)
		{
		return id2individuals.get(s);
		}
	@Override
	public java.util.Collection<? extends Sample> getSamples()
		{
		return this.id2individuals.values();
		}
	@Override
	public int hashCode() {
		return id.hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
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
	
	}




private static class PersonImpl implements Sample
		{
		private final FamilyImpl family;
		private final String id;
		String fatherId=null;
		String motherId=null;
		Sex sex=Sex.unknown;
		Status status=Status.unaffected;
		
		PersonImpl(final FamilyImpl family,final String id) {
			this.family = family;
			this.id = id;
			}
		
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
		
		private Sample getParent(final String s)
			{
			if(s==null || s.isEmpty() || s.equals("0")) return null;
			return getFamily().getSampleById(s);
			}
		
		public Sample getParent( int zeroOrOne) {
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
		
		public Sample getFather()
			{
			return getParent(fatherId);
			}
		
		public Sample getMother()
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
		public Sample validate() throws IllegalStateException {
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
}
