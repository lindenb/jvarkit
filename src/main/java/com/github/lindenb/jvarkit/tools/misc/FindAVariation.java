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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
/** 
 BEGIN_DOC
 
## Example

```
$ find ./ -name "*.vcf" -o -name "*.vcf.gz" |\
   java -jar dist/findavariation.jar -p "chr1:1234" 


htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12878	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12891	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12892	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12878	HOM_REF	C C
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12891	HET	C T
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12892	HET	C T
```

## History

  * 20180914 : replace DP4 with AD
 
 END_DOC
 */
@Program(name="findavariation",
	description="Finds a specific mutation in a list of VCF files",
	keywords={"vcf","variation","search"},
	modificationDate="20190409"
	)
public class FindAVariation extends Launcher
	{
	private static final Logger LOG = Logger.build(FindAVariation.class).make();
	@Parameter(names={"-p","--position"},description="A list of 'chrom/position'")
	private Set<String> positionsList = new HashSet<>();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-f","--posfile"},description="Add this file containing chrom:position")
	private Set<String> positionFilesList = new HashSet<>();
	@Parameter(names={"-homref","--homref"},description="Hide HOM_REF genotypes")
	private boolean hideHomRef=false;
	@Parameter(names={"-nocall","--nocall"},description="Hide NO_CALL genotypes")
	private boolean hideNoCall=false;
	@Parameter(names={"-snp","--snp"},description="Search only variant have the very same position (ignore overlapping variants)")
	private boolean onlySnp=false;
	@Parameter(names={"-indexed","--indexed"},description="[20171020] Search only in indexed vcf")
	private boolean indexedOnly=false;

	
	private static class Mutation
		{
		final String chrom;
		final int pos;
		Mutation(final String chrom,final int pos)
			{
			this.chrom=chrom;
			this.pos=pos;
			}
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + chrom.hashCode();
			result = prime * result + pos;
			return result;
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)return true;
			final Mutation other = (Mutation) obj;
			if (pos != other.pos) return false;
			 if (!chrom.equals(other.chrom))
				return false;
			
			return true;
		}
		
		@Override
		public String toString() {
			return chrom+":"+pos;
			}
		
		}
	private final Set<Mutation> mutations=new HashSet<Mutation>();
	private PrintWriter out=null;
	
    public FindAVariation()
    	{
    	}		
   
    
    private void reportPos(final String f,final VCFHeader header,final VariantContext ctx)
		{
		out.print(f);
		out.print('\t');
		out.print(ctx.getContig());
		out.print('\t');
		out.print(ctx.getStart());
		out.print('\t');
		out.print(ctx.getEnd());
		out.print('\t');
		out.print(ctx.hasID()?ctx.getID():".");
		out.print('\t');
		out.print(ctx.getReference().getDisplayString());
		}	

    
    private void report(
    		final String f,
    		final VCFHeader header,
    		final VariantContext ctx,
    		final Mutation mpos
    		)
    	{
    	
    	
    	if(!ctx.hasGenotypes())
    		{
    		reportPos(f,header,ctx);
    		out.println();
    		}
    	else
    		{
    		final GenotypesContext genotypes=ctx.getGenotypes();
    		for(int i=0;i< genotypes.size();++i)
    			{
    			Genotype g=genotypes.get(i);
    			if(!g.isCalled() && this.hideNoCall) continue;
    			if(g.isHomRef() && this.hideHomRef) continue;
    			reportPos(f,header,ctx);
    			out.print('\t');
    			out.print(g.getSampleName());
    			out.print('\t');
    			out.print(g.getType());
    			out.print('\t');
    			out.print(g.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(" ")));
    			out.print('\t');
    			if(g.hasAD())
    				{
					out.print(Arrays.stream(g.getAD()).mapToObj(x->String.valueOf(x)).collect(Collectors.joining(",")));//it's a String not an int[] ??
    				}
    			else
    				{
    				out.print('.');
    				}
    			out.print('\t');
    			out.print(g.hasDP()?String.valueOf(g.getDP()):".");
    			out.print('\t');
    			out.print(g.hasGQ()?String.valueOf(g.getGQ()):".");
    			out.println();
    			}
    		}
    	}	
    
    private Set<Mutation> convertFromVcfHeader(final String f,final VCFHeader h)
    	{
    	final SAMSequenceDictionary dict= h.getSequenceDictionary();
    	if(dict==null || dict.isEmpty()) {
    		return this.mutations;
    		}
    	
    	final Set<Mutation> copy=new HashSet<Mutation>(this.mutations.size());
    	final ContigNameConverter ctgCvt = ContigNameConverter.fromOneDictionary(dict);
    	
    	for(final Mutation m:this.mutations)
    		{
    		final String s=ctgCvt.apply(m.chrom);
    		if(StringUtils.isBlank(s))
    			{
    			LOG.warn("Cannot convert chrom "+m.chrom+" in "+f);
    			continue;
    			}
    		copy.add(new Mutation(s, m.pos));
    		}
    	return copy;
    	}

    
    
    private void scanLine(final String line)  {
    	if(StringUtils.isBlank(line) || line.startsWith("#")) return;
    	if(IOUtil.isUrl(line)) {
    		scanRemote(line);
    		}
    	else
    		{
    		scanPath(line);
    		}
    	}
    private void scanRemote(final String url)  {
    	TabixVcfFileReader tabix = null;
    	try {
    		tabix=new TabixVcfFileReader(url);
    		final VCFHeader header = tabix.getHeader();
    		final Set<String> chromosomes = tabix.getChromosomes();
    		for(final Mutation m:convertFromVcfHeader(url,header))
    			{
    			if(!chromosomes.contains(m.chrom)) continue;
    			final java.util.Iterator<VariantContext> iter2 = tabix.iterator(m.chrom, m.pos, m.pos);
    		    while(iter2.hasNext())
                      {
                      final VariantContext ctx=iter2.next();
                      if(this.onlySnp )
                            {       
                            if(ctx.getStart()!=m.pos || ctx.getEnd()!=m.pos) continue;
                            }
                      report(url,header,ctx,m);
                      }
                  CloserUtil.close(iter2);
                  }
    		tabix.close();
    		tabix=null;
    		}
    	catch(final htsjdk.tribble.TribbleException.InvalidHeader err)
			{
			LOG.warn(url+"\t"+err.getMessage());
			}
		catch(final Throwable err)
			{
			LOG.severe("cannot read "+url,err);
			}
    	finally
    		{
    		CloserUtil.close(tabix);
    		}
    	}

    
    
    private void scanPath(final String vcfPathString)  {
    	final Path vcfPath=Paths.get(vcfPathString);
    	if(Files.isDirectory(vcfPath)) return;
    	if(!Files.isReadable(vcfPath)) return;
		VCFIterator iter=null;
		VCFFileReader r=null;
		try {
			if(VCFUtils.isTribbleVcfPath(vcfPath) || VCFUtils.isTribbleVcfPath(vcfPath))
				{
				r=new VCFFileReader(vcfPath,true);
				final VCFHeader header =r.getFileHeader();
				for(final Mutation m:convertFromVcfHeader(vcfPath.toString(),header))
					{
					try( CloseableIterator<VariantContext> iter2 = r.query(m.chrom, m.pos, m.pos)) {
						while(iter2.hasNext())
							{
							final VariantContext ctx=iter2.next();
							if(this.onlySnp )
								{	
								if(ctx.getStart()!=m.pos || ctx.getEnd()!=m.pos) continue;
								}
							report(vcfPathString,header,ctx,m);
							}
						}
					}
				r.close();
				r=null;
				}   				
			else if(!this.indexedOnly)
				{
				iter=VCFUtils.createVCFIteratorFromPath(vcfPath);
				final VCFHeader header = iter.getHeader();
				final Set<Mutation> mutlist=convertFromVcfHeader(vcfPath.toString(),iter.getHeader());
				while(iter.hasNext())
					{
					final VariantContext ctx=iter.next();
					final Mutation m=new Mutation(ctx.getContig(), ctx.getStart());
					
					for(final Mutation m2: mutlist)
						{
						if(m.equals(m2)) {
					    	if(this.onlySnp )
								{	
								if(ctx.getStart()!=m2.pos || ctx.getEnd()!=m2.pos) continue;
								}	
							report(vcfPathString,header,ctx,m2);
							break;
							}
						}
					}
				iter.close();
				iter=null;
				}
			}
		catch(final htsjdk.tribble.TribbleException.InvalidHeader err)
			{
			LOG.warn(vcfPathString+"\t"+err.getMessage());
			}
		catch(final Throwable err)
			{
			LOG.severe("cannot read "+vcfPathString,err);
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(iter);
			}
	}
    

	private Mutation parseMutation(final String s)
		{
		final int colon=s.indexOf(':');
		if(colon==-1 || colon+1==s.length())
			{
			throw new IllegalArgumentException("Bad chrom:pos "+s);
			}
		
		final String chrom=s.substring(0,colon).trim();
		if(StringUtils.isBlank(chrom))
			{
			throw new IllegalArgumentException("Bad chrom:pos "+s);
			}
		final Mutation m=new Mutation(chrom, Integer.parseInt(s.substring(colon+1)));
		return m;
		}
	
	@Override
	public int doWork(final List<String> args) {
		BufferedReader r=null;
		try
			{
			for(final String f:this.positionFilesList)
				{
				r = IOUtils.openURIForBufferedReading(f);
				String line;
				while((line=r.readLine())!=null)
					{
					if(line.isEmpty() || line.startsWith("#")) continue;
					final Mutation m= parseMutation(line);
					LOG.debug("adding "+m);
					this.mutations.add(m);
					}
				r.close();
				}
			
			for(final String s:this.positionsList)
				{
				final Mutation m= parseMutation(s);
				LOG.debug("adding "+m);
				this.mutations.add(m);
				}			
			
			this.out=super.openFileOrStdoutAsPrintWriter(this.outputFile);
			this.out.println("#FILE\tCHROM\tstart\tend\tID\tREF\tsample\ttype\tALLELES\tAD\tDP\tGQ");
			
			
			
			if(args.isEmpty())
				{
				BufferedReader br= new BufferedReader(new InputStreamReader(stdin()));
				br.lines().forEach(L->scanLine(L));
				br.close();
				}
			else if(args.size()==1 && args.get(0).endsWith(".list"))
				{
				BufferedReader br= IOUtils.openPathForBufferedReading(Paths.get(args.get(0)));
				br.lines().forEach(L->scanLine(L));
				br.close();
				}
			else 
				{
				args.forEach(L->scanLine(L));
				}
			this.out.flush();
			this.out.close();
			this.out=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(r);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(final String[] args) {
		new FindAVariation().instanceMainWithExit(args);

	}

}
