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
package com.github.lindenb.jvarkit.pedigree;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.File;
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



public class PedigreeParser {

/** default sex parser for string */ 
private static final Function<String, Sex> DEFAULT_SEX_PARSER = S->{
	if(S.length()==1) {
		switch(S.charAt(0)) {
			case '0': return Sex.unknown;
			case '1': return Sex.male;
			case '2': return Sex.female;
			}
		}
	throw new IllegalArgumentException("bad sex in \""+S+"\", should be 0/1/2.");
	};

/* parse sex from string */
private Function<String, Sex> sexParser = DEFAULT_SEX_PARSER;

public PedigreeParser setSexParser(final Function<String, Sex> fun) {
	this.sexParser = fun;
	return this;
	}

public Function<String,Sex>  getSexParser() {
	return this.sexParser;
	}

/** default status parser for string */ 
private static final Function<String,Status> DEFAULT_STATUS_PARSER = S->{
	if(S.equals("-9")) return Status.missing;
	if(S.equals("0")) return Status.unaffected;
	if(S.equals("1")) return Status.affected;
	throw new IllegalArgumentException("bad status in \""+S+"\", should be 0/1/-9 .");
	};

/* parse status from string */
private Function<String,Status> statusParser = DEFAULT_STATUS_PARSER;

public PedigreeParser setStatusParser(final Function<String,Status> fun) {
	this.statusParser = fun;
	return this;
	}

public Function<String,Status>  getStatusParser() {
	return this.statusParser ;
	}

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
	@Override
	public String toString() {
		return "pedigree families.count=" + getFamilies().size();
		}
	}

/** parse pedigree file */
public Pedigree parse(final File pedFile) throws IOException {
	try(final BufferedReader br=IOUtils.open(pedFile)) {
		return parse(br);
		}
	}

/** parse pedigree file */
public Pedigree parse(final BufferedReader br)throws IOException
	{
	final CharSplitter tab = CharSplitter.TAB;
	final PedigreeImpl ped = new PedigreeImpl();
	String line;
	while((line=br.readLine())!=null) {
		if(StringUtils.isBlank(line)) continue;
		if(line.startsWith("#")) continue;
		final String tokens[]= tab.split(line);
		if(tokens.length<4) throw new IllegalArgumentException("not enough tokens for pedigree in "+line);
		final String famId= tokens[0];
		final String indiId = tokens[1];
		final String fatherId = tokens[2];
		final String motherId = tokens[3];
		final String sex = (tokens.length>4?tokens[4]:"");
		final String status = (tokens.length>5?tokens[5]:"");
		build(ped, famId, indiId, fatherId, motherId, sex, status);
		}
	return ped.validate();
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
	
	if(!StringUtils.isBlank(sexxx))
		{
		p.sex = getSexParser().apply(sexxx);
		}
	
	if(!StringUtils.isBlank(status))
		{
		p.status = getStatusParse().apply(status);
		}
	else if(this.statusRequired) {
		throw new IllegalArgumentException("status must be declared");
		}
		
	fam.id2individuals.put(p.id, p);		
	}

private static class FamilyImpl implements Family
	{
	private final PedigreeImpl ped;
	private final String id;
	private Map<String,SampleImpl> id2individuals = new TreeMap<>();
	private FamilyImpl(final PedigreeImpl ped,final String id) {
		this.ped = ped;
		this.id = id;
		if(StringUtils.isBlank(id)) throw new IllegalArgumentException("Family id cannot be empty");
		}

	@Override
	public String getId() {
		return this.id;
		}
	@Override
	public Sample getSampleById(final String s) {
		return this.id2individuals.get(s);
		}
	@Override
	public java.util.Collection<? extends Sample> getSamples()
		{
		return Collections.unmodifiableCollection(this.id2individuals.values());
		}
	@Override
	public int hashCode() {
		return id.hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if (this == obj) return true;
		if (obj == null || !(obj instanceof FamilyImpl)) {
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
			if(StringUtils.isBlank(id)) throw new IllegalArgumentException("bad sample empty id");
			if(id.equals("0")) throw new IllegalArgumentException("sample cannot be named id");
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
			} }
		
		private boolean hasParent(final String s)
			{
			return !(StringUtils.isBlank(s) || s.equals("0"));
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
		public boolean equals(final Object obj) {
			if (this == obj) return true;
			if (obj == null || !(obj instanceof SampleImpl)) {
				return false;
			}
			final PersonImpl p =(PersonImpl)obj;
			return	this.family.equals(p.family) &&
					this.id.equals(p.id);
			}
		
		//@Override
		public boolean hasUniqId() 
			{
			return (Pedigree.this.families.values().
					stream().
					flatMap(F->F.getIndividuals().stream())).
					filter(P->getId().equals(P.getId())).count() ==1L;
			}
		
		//@Override
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
