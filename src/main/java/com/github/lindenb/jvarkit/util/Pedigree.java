/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
import java.util.HashSet;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;

import com.github.lindenb.jvarkit.io.IOUtils;


public class Pedigree
	{
	private Map<String,FamilyImpl > families=new TreeMap<String, Pedigree.FamilyImpl>();
	
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
	
	public interface Family
		{
		public String getId();
		public Person getPersonById(String s);
		public java.util.Collection<? extends Person> getIndividuals();
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
		}
	
	public interface Person
		{
		public String getId();
		public Family getFamily();
		public Person getFather();
		public Person getMother();
		public Sex getSex();
		public Status getStatus();
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
		
		private Person getParent(String s)
			{
			if(s==null || s.isEmpty() || s.equals("0")) return null;
			return getFamily().getPersonById(s);
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
		public String toString() {
			return family+":"+this.id;
			}
		}
	
	private Pedigree()
		{
		
		}
	/** get all the families in this pedigree */
	public java.util.Collection<? extends Family> getFamilies()
		{
		return this.families.values();
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
	
	private void read(final BufferedReader r) throws IOException
		{
		Pattern tab = Pattern.compile("[\t]");
		String line;
		while((line=r.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			final String tokens[]=tab.split(line);
			FamilyImpl fam=this.families.get(tokens[0]);
			if(fam==null)
				{
				fam=new FamilyImpl();
				this.families.put(tokens[0], fam);
				}
			if(fam.getPersonById(tokens[1])!=null) throw new IOException("duplicate individual: "+line);
			PersonImpl p=new PersonImpl();
			p.family=fam;
			p.id=tokens[1];
			p.fatherId=tokens[2];
			p.motherId=tokens[3];
			
			if(tokens.length>4)
				{
				if(tokens[4].equals("1")) p.sex=Sex.male;
				else if(tokens[4].equals("2")) p.sex=Sex.female;
				}
			
			if(tokens.length>5)
				{
				if(tokens[5].equals("1")) p.status=Status.affected;
				else if(tokens[5].equals("0")) p.status=Status.unaffected;
				}
			
			fam.individuals.put(p.id, p);
			}
		}
	
	public static Pedigree readPedigree(final File f) throws IOException
		{
		final BufferedReader r= IOUtils.openFileForBufferedReading(f);
		final Pedigree ped=Pedigree.readPedigree(r);
		r.close();
		return ped;
		}	
	public static Pedigree readPedigree(final BufferedReader r) throws IOException
		{
		final Pedigree ped=new Pedigree();
		ped.read(r);
		return ped;
		}	
	}
