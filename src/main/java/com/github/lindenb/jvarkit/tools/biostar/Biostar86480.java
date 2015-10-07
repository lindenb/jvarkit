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

import htsjdk.samtools.util.CloserUtil;

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;



import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.command.Command;

public class Biostar86480 extends AbstractBiostar86480
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Biostar86480.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractBiostar86480.AbstractBiostar86480Command
		{    

		private Rebase rebase=Rebase.createDefaultRebase();
	
	
	
	private void digest(
			String seqName,
			int position0,
			final List<Character> sequence,
			PrintStream out
			)
		{
		for(Rebase.Enzyme enzyme:this.rebase)
			{
			if(enzyme.size()>sequence.size()) continue;
			for(int strand=0;strand<2;++strand)
				{
				int x=0;
				for(x=0;x< enzyme.size();++x)
					{
					char c=(strand==0?
							enzyme.at(x):
							AcidNucleics.complement(enzyme.at((enzyme.size()-1)-x))
							);
					if(!Rebase.compatible(sequence.get(x),c)) break;
					}
				if(x==enzyme.size())
					{
					out.print(seqName);
					out.print('\t');
					out.print(position0);
					out.print('\t');
					out.print(position0+enzyme.size());
					out.print('\t');
					for(int y=0;y< enzyme.size();++y)
						{
						out.print(sequence.get(y));
						}
					out.print('\t');
					out.print(1000);
					out.print('\t');
					out.print(strand==1?'-':'+');
					out.print('\t');
					out.print(enzyme.getName());
					out.print('\t');
					out.print(enzyme.getDecl());
					out.println();
					break;
					}
				if(enzyme.isPalindromic()) break;
				}
			}
		}
	
	private void run(Reader in,PrintStream out) throws IOException
		{
		int longest=0;
		for(Rebase.Enzyme E:this.rebase)
			{
			longest=Math.max(E.size(), longest);
			}
		String seqName="";
		int position0=0;
		ArrayList<Character> sequences=new ArrayList<Character>(longest);
		for(;;)
			{
			int c=in.read();
			if(c==-1 || c=='>')
				{
				while(!sequences.isEmpty())
					{
					digest(seqName,position0,sequences,out);
					++position0;
					sequences.remove(0);
					}
				if(c==-1) break;
				StringBuilder b=new StringBuilder();
				while((c=in.read())!=-1 && c!='\n')
					{
					b.append((char)c);
					}
				seqName=b.toString();
				position0=0;
				}
			else if(!Character.isWhitespace(c))
				{
				sequences.add((char)Character.toUpperCase(c));
				if(sequences.size()==longest)
					{
					digest(seqName,position0,sequences,out);
					++position0;
					sequences.remove(0);
					if(position0%1000000==0)
						{
						LOG.info(seqName+" "+position0);
						}
					}
				}
			}
		}

	@Override
	public Collection<Throwable> call() throws Exception
		{
		final List<String> args = getInputFiles();
		if(!super.onlyEnz.isEmpty())
			{
			Rebase rebase2=new Rebase();
			for(String e:super.onlyEnz)
				{
				Rebase.Enzyme enz=this.rebase.getEnzymeByName(e);
				if(enz==null)
					{
					LOG.error("Cannot find enzyme "+enz +" in RE list.");
					LOG.error("Current list is:");
					for(Rebase.Enzyme E: this.rebase)
						{
						LOG.error("\t"+E);
						}
					return wrapException("Cannot find enzyme "+enz +" in RE list.");
					}
				rebase2.getEnzymes().add(enz);
				}
			this.rebase=rebase2;
			}
		Reader in=null;
		PrintStream out=null;
		try
			{
			if(getOutputFile()==null)
				{
				out=stdout();
				}
			else
				{
				out=new PrintStream(getOutputFile());
				}
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				in = new InputStreamReader(stdin());
				run(in,out);
				in.close();
				}
			else
				{
				for(String arg: args)
					{
					LOG.info("Opening "+arg);
					in=IOUtils.openURIForBufferedReading(arg);
					run(in,out);
					in.close();
					}
				}
			return Collections.emptyList();
			}
		catch(Throwable err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(out);
			}
		}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new Biostar86480().instanceMainWithExit(args);
		}

	}
