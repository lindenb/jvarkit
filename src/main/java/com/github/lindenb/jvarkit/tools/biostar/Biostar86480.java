/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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

import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.Path;
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
	keywords={"rebase","genome","enzyme","restricion","genome"},
	creationDate="20131114",
	modificationDate="20210805"
	)
public class Biostar86480 extends Launcher
	{
	private static final Logger LOG = Logger.build(Biostar86480.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-E","--enzyme"},description="restrict to that enzyme name.")
	private Set<String> onlyEnz = new HashSet<>();
	@Parameter(names={"--min-size","--min-weight"},description="restrict to that enzyme 'size/weight'. ignore if 'x' <=0")
	private float min_size= 0;
	@Parameter(names={"-l"},description="list available enzymes",help=true)
	private boolean dump_enzymes = false;


	private final Rebase rebase=Rebase.createDefaultRebase();
	
	public Biostar86480()
		{
		}

	
	
	private void digest(
			final String seqName,
			final int position0,
			final CharSequence sequence,
			final PrintWriter out
			)
		{
		for(final Rebase.Enzyme enzyme:this.rebase)
			{
			if(enzyme.size()>sequence.length()) continue;
			for(int strand=0;strand<2;++strand)
				{
				int x=0;
				for(x=0;x< enzyme.size();++x)
					{
					char c=(strand==0?
							enzyme.at(x):
							AcidNucleics.complement(enzyme.at((enzyme.size()-1)-x))
							);
					if(!Rebase.compatible(sequence.charAt(x),c)) break;
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
						out.print(sequence.charAt(y));
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
	
	private void run(final Reader in,final PrintWriter out) throws IOException
		{
		final int longest= this.rebase.stream().mapToInt(E->E.size()).max().orElse(0);
		
		String seqName="";
		int position0=0;
		final StringBuilder sequences=new StringBuilder(longest);
		for(;;)
			{
			int c=in.read();
			if(c==-1 || c=='>')
				{
				while(sequences.length()>0)
					{
					digest(seqName,position0,sequences,out);
					++position0;
					sequences.delete(0, 1);
					}
				if(c==-1) break;
				final StringBuilder b=new StringBuilder();
				while((c=in.read())!=-1 && c!='\n')
					{
					b.append((char)c);
					}
				seqName=b.toString();
				position0=0;
				}
			else if(!Character.isWhitespace(c))
				{
				sequences.append((char)Character.toUpperCase(c));
				if(sequences.length()==longest)
					{
					digest(seqName,position0,sequences,out);
					++position0;
					sequences.delete(0, 1);
					}
				}
			}
		}

	@Override
	public int doWork(final List<String> args) {
		
		
		if(!onlyEnz.isEmpty())
			{
			for(String e:this.onlyEnz)
				{
				final Rebase.Enzyme enz=this.rebase.getEnzymeByName(e);
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
				}
			this.rebase.removeIf(ENZ->onlyEnz.contains(ENZ.getName()));
			}
		
		if(this.min_size>0f) {
			this.rebase.removeIf(ENZ->ENZ.getWeight()< this.min_size);
			}
		
		if(dump_enzymes) {
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				out.println("#NAME,SITE,SIZE,WEIGHT");
				for(final Rebase.Enzyme enz:this.rebase) {
					out.println(enz.getName()+","+enz.getDecl()+","+enz.size()+","+enz.getWeight());
				}
				out.flush();
				return 0;
			}
			catch(IOException err) {
				LOG.error(err);
				return -1;
			}
		}

		
		if(this.rebase.isEmpty()) {
			LOG.error("No enzyme in database");
			return -1;
			}
		try
			{
			try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				if(args.isEmpty())
					{
					try(Reader in=new InputStreamReader(stdin())) {
						run(in,out);
						}
					}
				else
					{
					for(final String arg:args)
						{
						try(Reader in=IOUtils.openURIForBufferedReading(arg)) {
							run(in,out);
							}
						}
						
					}
				out.flush();
				}
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
