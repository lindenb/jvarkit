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


*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;

import htsjdk.samtools.util.CloserUtil;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.command.Command;
import com.sleepycat.bind.tuple.StringBinding;
import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;
import com.sleepycat.je.Cursor;
import com.sleepycat.je.Database;
import com.sleepycat.je.DatabaseConfig;
import com.sleepycat.je.DatabaseEntry;
import com.sleepycat.je.Environment;
import com.sleepycat.je.EnvironmentConfig;
import com.sleepycat.je.LockMode;
import com.sleepycat.je.OperationStatus;
import com.sleepycat.je.Transaction;

/**
 * @author lindenb
 *
 */
public class Biostar92368 extends AbstractBiostar92368
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar95652.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar92368.AbstractBiostar92368Command
		{    

	
	private Environment environment=null;
	private Database database=null;
	
	private class StringSetBinding extends TupleBinding<Set<String>>
		{
		@Override
		public Set<String> entryToObject(TupleInput in)
			{
			int n=in.readInt();
			Set<String> h=new HashSet<String>(n);
			for(int i=0;i< n;++i) h.add(in.readString());
			return h;
			}
		public void objectToEntry(Set<String> h, TupleOutput out)
			{
			out.writeInt(h.size());
			for(String s:h) out.writeString(s);
			}
		}
	private StringSetBinding bindingSet=new StringSetBinding();
	
	
	private int recursive(
			Transaction txn,
			final String endProt,
			Stack<String> path,
			int bestDepth,
			final int maxDepth
			)
		{
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		StringBinding.stringToEntry(path.lastElement(), key);
		if(this.database.get(txn,key,data,null)==OperationStatus.SUCCESS)
			{
			Set<String> partners=this.bindingSet.entryToObject(data);
			for(String prot2:partners)
				{
				int depth=-1;
				if(prot2.equals(endProt))
					{
					depth=path.size()-1;
					}
				else if(path.size()<maxDepth && !path.contains(prot2))
					{
					path.push(prot2);
					depth=recursive(txn,endProt,path,bestDepth,maxDepth);
					path.pop();
					}
				
				if(depth!=-1 && (bestDepth==-1 || bestDepth>depth))
					{
					bestDepth=depth;
					}
				}
			}
		return bestDepth;
		}
	
	private void load(Transaction txn,LineIterator r)
		{
		long nLines=0L;
		DatabaseEntry key=new DatabaseEntry();
		DatabaseEntry data=new DatabaseEntry();
		Pattern tab=Pattern.compile("[\t ]");
		Cursor c=this.database.openCursor(txn, null);
		while(r.hasNext())
			{
			if(nLines++ % 10000 ==0)
				{
				LOG.info("Lines : "+nLines);
				}
			String line=r.next();
			if(line.isEmpty() || line.startsWith("#")) continue;
			String tokens[]=tab.split(line);
			for(int i=0;i< 2;++i)
				{
				String a=tokens[i==0?0:1];
				String b=tokens[i==0?1:0];
				StringBinding.stringToEntry(a, key);
				OperationStatus status=c.getSearchKey(key, data, LockMode.DEFAULT);
				Set<String> partners; 
				if(status==OperationStatus.SUCCESS)
					{
					partners=bindingSet.entryToObject(data);
					
					}
				else
					{
					partners=new HashSet<String>(1);
					}
				partners.add(b);
				bindingSet.objectToEntry(partners,data);
				if(status==OperationStatus.SUCCESS)
					{
					status=c.putCurrent(data);
					}
				else
					{
					status=c.put(key, data);
					}
				if(status!=OperationStatus.SUCCESS)
					{
					throw new RuntimeException("Cannot insert data in bdb "+a+"/"+b);
					}
				if(a.equals(b)) break;//self-self
				}
			}
		CloserUtil.close(c);
		}
	
	

	
	
		@Override
		public Collection<Throwable> call() throws Exception
			{
			
			if(dbHome==null)
				{
				return wrapException("Undefined DB-Home");
				}
			environment=null;
			EnvironmentConfig envCfg=new EnvironmentConfig();
			this.database=null;
			Transaction txn=null;
			try
				{
				envCfg.setAllowCreate(true);
				this.environment=new Environment(dbHome, envCfg);
				DatabaseConfig cfg=new DatabaseConfig();
				cfg.setAllowCreate(true);
				cfg.setTemporary(true);
				this.database=environment.openDatabase(txn, "interactions", cfg);
				List<String> args=getInputFiles();
				if(args.isEmpty())
					{
					LOG.info("reading stdin");
					LineIterator r=IOUtils.openStdinForLineIterator();
					load(txn,r);
					CloserUtil.close(r);
					}
				else
					{
					for(String filename : args)
						{
						LOG.info("reading "+filename);
						LineIterator r=IOUtils.openURIForLineIterator(filename);
						load(txn,r);
						CloserUtil.close(r);
						}
					}
				DatabaseEntry key1=new DatabaseEntry();
				DatabaseEntry data1=new DatabaseEntry();
				DatabaseEntry key2=new DatabaseEntry();
				DatabaseEntry data2=new DatabaseEntry();
	
				Cursor c1=this.database.openCursor(txn, null);
				Cursor c2=this.database.openCursor(txn, null);
				while(c1.getNext(key1, data1, LockMode.DEFAULT)==OperationStatus.SUCCESS)
					{
					
					String prot1=StringBinding.entryToString(key1);
					
					boolean first=true;
					
					while((first?c2.getFirst(key2, data2, LockMode.DEFAULT):c2.getNext(key2, data2, LockMode.DEFAULT))==OperationStatus.SUCCESS)
						{
						first=false;
						String prot2=StringBinding.entryToString(key2);
						if(prot2.compareTo(prot1)<=0) continue;
						
						Stack<String> path=new Stack<String>();
						path.push(prot1);
						int depth=recursive(txn,prot2,path,-1,maxDepth);
						if(depth!=-1)
							{
							System.out.println(prot1+"\t"+prot2+"\t"+depth);
							}
						else
							{
							//System.out.println(prot1+"\t"+prot2+"\t"+depth);
							}
						if(System.out.checkError()) break;
						}
					
					}
				CloserUtil.close(c2);
				CloserUtil.close(c1);
				return Collections.emptyList();
				}
			catch (Exception err)
				{
				LOG.error(err);
				return wrapException(err);
				}
			finally
				{
				CloserUtil.close(this.database);
				CloserUtil.close(this.environment);
				}
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar92368().instanceMainWithExit(args);

		}

	}
