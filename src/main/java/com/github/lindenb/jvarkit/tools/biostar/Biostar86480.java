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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.biostar;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.Rebase;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**
BEGIN_DOC

## Example
```bash
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz" |\
gunzip -c  |\
java -jar dist/biostar86480.jar -E AarI -E EcoRI  

chr3	60645	60651	GAATTC	1000	+	EcoRI	G^AATTC
chr3	60953	60959	GAATTC	1000	+	EcoRI	G^AATTC
chr3	68165	68172	GCAGGTG	1000	-	AarI	CACCTGC(4/8)
chr3	70263	70269	GAATTC	1000	+	EcoRI	G^AATTC
chr3	70945	70952	GCAGGTG	1000	-	AarI	CACCTGC(4/8)
chr3	71140	71146	GAATTC	1000	+	EcoRI	G^AATTC
chr3	72264	72270	GAATTC	1000	+	EcoRI	G^AATTC
chr3	74150	74156	GAATTC	1000	+	EcoRI	G^AATTC
chr3	75063	75069	GAATTC	1000	+	EcoRI	G^AATTC
chr3	78438	78444	GAATTC	1000	+	EcoRI	G^AATTC
chr3	81052	81059	CACCTGC	1000	+	AarI	CACCTGC(4/8)
chr3	84498	84504	GAATTC	1000	+	EcoRI	G^AATTC
chr3	84546	84552	GAATTC	1000	+	EcoRI	G^AATTC
chr3	84780	84787	CACCTGC	1000	+	AarI	CACCTGC(4/8)
chr3	87771	87777	GAATTC	1000	+	EcoRI	G^AATTC
chr3	95344	95351	GCAGGTG	1000	-	AarI	CACCTGC(4/8)
chr3	96358	96364	GAATTC	1000	+	EcoRI	G^AATTC
chr3	96734	96740	GAATTC	1000	+	EcoRI	G^AATTC
chr3	105956	105962	GAATTC	1000	+	EcoRI	G^AATTC
chr3	107451	107457	GAATTC	1000	+	EcoRI	G^AATTC
(...)
```
END_DOC
 */
@Program(name="biostar86480",description="Genomic restriction finder",
	biostars=86480,
		keywords={"rebase","genome","enzyme","restricion","genome"}
		)
public class Biostar86480 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar86480.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-E","--enzyme"},description="restrict to that enzyme.")
	private Set<String> onlyEnz = new HashSet<>();

	private Rebase rebase=Rebase.createDefaultRebase();
	
	public Biostar86480()
		{
		}

	
	
	private void digest(
			final String seqName,
			int position0,
			final List<Character> sequence,
			final PrintStream out
			)
		{
		for(final Rebase.Enzyme enzyme:this.rebase)
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
	public int doWork(final List<String> args) {
		
		if(!onlyEnz.isEmpty())
			{
			Rebase rebase2=new Rebase();
			for(String e:onlyEnz)
				{
				Rebase.Enzyme enz=this.rebase.getEnzymeByName(e);
				if(enz==null)
					{
					LOG.error("Cannot find enzyme "+e +" in RE list.");
					LOG.error("Current list is:");
					for(final Rebase.Enzyme E: this.rebase)
						{
						LOG.error("\t"+E);
						}
					return -1;
					}
				rebase2.getEnzymes().add(enz);
				}
			this.rebase=rebase2;
			}
		PrintStream out;
		try
			{
			out = super.openFileOrStdoutAsPrintStream(this.outputFile);
			if(args.isEmpty())
				{
				LOG.info("Reading from stdin");
				run(new InputStreamReader(stdin()),out);
				}
			else
				{
				for(final String arg:args)
					{
					LOG.info("Opening "+arg);
					final Reader in=IOUtils.openURIForBufferedReading(arg);
					run(in,out);
					in.close();
					}
					
				}
			out.flush();
			out.close();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
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
