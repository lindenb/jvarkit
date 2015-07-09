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
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;


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
		}
	
	private Pedigree()
		{
		
		}
	
	public java.util.Collection<? extends Family> getFamilies()
		{
		return this.families.values();
		}
	
	private void read(BufferedReader r) throws IOException
		{
		String line;
		while((line=r.readLine())!=null)
			{
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=line.split("\t");
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
	
	public static Pedigree readPedigree(File f) throws IOException
		{
		BufferedReader r=new BufferedReader(new FileReader(f));
		Pedigree ped=Pedigree.readPedigree(r);
		r.close();
		return ped;
		}	
	public static Pedigree readPedigree(BufferedReader r) throws IOException
		{
		Pedigree ped=new Pedigree();
		ped.read(r);
		return ped;
		}	
	}
