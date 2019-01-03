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
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.regex.Pattern;

import htsjdk.tribble.readers.LineIterator;

import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
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
BEGIN_DOC

##Example

```bash
$ cat input.txt
A1	A2
A2	A3
A3	A4
A4	A5
A1	A6


$ mkdir -p tmp
$  java -jar dist/biostar92368.jar -D tmp input.txt

A1	A2	0
A1	A3	1
A1	A4	2
A1	A6	0
A2	A3	0
A2	A4	1
A2	A5	2
A2	A6	1
A3	A4	0
A3	A5	1
A3	A6	2
A4	A5	0

```
END_DOC

 *
 */
@Program(name="biostar92368",
	biostars=92368,
	keywords={"protein","interaction","interactome"},
	description="Binary interactions depth.")
public class Biostar92368 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar92368.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-D","--bdbhome"},description="berkeleydb home",required=true)
	private File dbHome=null;
	@Parameter(names={"-M","--maxdepth"},description="Max depth")
	private int maxDepth=3;

	
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
	public int doWork(List<String> args) {
		if(this.dbHome==null)
			{
			LOG.error("Undefined DB-Home");
			return -1;
			}
		environment=null;
		EnvironmentConfig envCfg=new EnvironmentConfig();
		this.database=null;
		Transaction txn=null;
		PrintStream out=null;
		try
			{
			
			envCfg.setAllowCreate(true);
			this.environment=new Environment(dbHome, envCfg);
			DatabaseConfig cfg=new DatabaseConfig();
			cfg.setAllowCreate(true);
			cfg.setTemporary(true);
			this.database=environment.openDatabase(txn, "interactions", cfg);
			if(args.isEmpty())
				{
				LOG.info("reading stdin");
				LineIterator r=IOUtils.openStdinForLineIterator();
				load(txn,r);
				CloserUtil.close(r);
				}
			else
				{
				for(String filename:args)
					{
					LOG.info("reading "+filename);
					LineIterator r=IOUtils.openURIForLineIterator(filename);
					load(txn,r);
					CloserUtil.close(r);
					}
				}
			out= super.openFileOrStdoutAsPrintStream(outputFile);
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
						out.println(prot1+"\t"+prot2+"\t"+depth);
						}
					else
						{
						//System.out.println(prot1+"\t"+prot2+"\t"+depth);
						}
					if(out.checkError()) break;
					}
				
				}
			CloserUtil.close(c2);
			CloserUtil.close(c1);
			out.flush();
			out.close();
			return 0;
			}
		catch (Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.database);
			CloserUtil.close(this.environment);
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
