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


*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
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
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.vcf.BcfToolsUtils;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;
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
	keywords={"vcf","variation","search","find","bcf"},
	creationDate="20140623",
	modificationDate="20200217"
	)
public class FindAVariation extends Launcher
	{
	private static final Logger LOG = Logger.build(FindAVariation.class).make();
	@Parameter(names={"-p","--position"},description="A list of 'chrom/position'")
	private Set<String> positionsList = new HashSet<>();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
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
	@Parameter(names={"--bcf"},description="Enable BCF. On 20200217 the java library htsjdk doesn't support the recent version of bcf. This tool will call `bcftools` to get the variants. bcftools must be in the ${PATH}.")
	private boolean use_bcf=false;

	
	private final Set<SimplePosition> mutations=new HashSet<>();
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
    		final SimplePosition mpos
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
    			final Genotype g=genotypes.get(i);
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
    
    private Set<SimplePosition> convertFromVcfHeader(final String f,final VCFHeader h)
    	{
    	final SAMSequenceDictionary dict= h.getSequenceDictionary();
    	if(dict==null || dict.isEmpty()) {
    		return this.mutations;
    		}
    	
    	final Set<SimplePosition> copy=new HashSet<>(this.mutations.size());
    	final ContigNameConverter ctgCvt = ContigNameConverter.fromOneDictionary(dict);
    	
    	for(final SimplePosition m:this.mutations)
    		{
    		final String s=ctgCvt.apply(m.getContig());
    		if(StringUtils.isBlank(s))
    			{
    			LOG.warn("Cannot convert chrom "+m.getContig()+" in "+f);
    			continue;
    			}
    		copy.add(m.renameContig(s));
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
    	VCFReader tabix = null;
    	try {
    		tabix=VCFReaderFactory.makeDefault().open(url);
    		final VCFHeader header = tabix.getHeader();
    		final SAMSequenceDictionary dict = header.getSequenceDictionary();
    		for(final SimplePosition m:convertFromVcfHeader(url,header))
    			{
    			if(dict.getSequence(m.getContig())==null) continue;
    			final java.util.Iterator<VariantContext> iter2 = tabix.query(new Interval(m));
    		    while(iter2.hasNext())
                      {
                      final VariantContext ctx=iter2.next();
				      if(this.onlySnp  && (ctx.getStart()!=m.getPosition() || ctx.getEnd()!=m.getPosition())) continue;
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
    private boolean isVcf(final Path f) {
    	final String filename=f.getFileName().toString();
    	return FileExtensions.VCF_LIST.stream().anyMatch(EXT->filename.endsWith(EXT)) ||
    		filename.endsWith(".vcf.bgz");
    	}
    private boolean isIndexed(final Path f) {
    	if(!isVcf(f)) return false;
    	final String filename=f.getFileName().toString();

    	for(final String suff:new String[]{FileExtensions.TRIBBLE_INDEX,FileExtensions.TABIX_INDEX,FileExtensions.CSI}) {
    	
			final Path index=Paths.get(
					f.getParent().toString(),
					filename+ suff
					);
			if( Files.exists(index) && !Files.isDirectory(index)) return true;
	    	}
    	return false;
    	}
    
    private void scanPath(final String vcfPathString)  {
    	final Path vcfPath=Paths.get(vcfPathString);
    	if(Files.isDirectory(vcfPath)) return;
    	if(!Files.isReadable(vcfPath)) return;
		VCFIterator iter=null;
		VCFReader r=null;
		try {
			if(vcfPathString.endsWith(FileExtensions.BCF) && BcfToolsUtils.isBcfToolsRequired(Paths.get(vcfPathString))) {
				if(!this.use_bcf ) return;
				}
			
			if(isIndexed(vcfPath))
				{
				r= VCFReaderFactory.makeDefault().open(vcfPath,true);
				final VCFHeader header =r.getHeader();
				for(final SimplePosition m:convertFromVcfHeader(vcfPath.toString(),header))
					{
					try( CloseableIterator<VariantContext> iter2 = r.query(m)) {
						while(iter2.hasNext())
							{
							final VariantContext ctx=iter2.next();
					    	if(this.onlySnp  && (ctx.getStart()!=m.getPosition() || ctx.getEnd()!=m.getPosition())) continue;
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
				final Set<SimplePosition> mutlist=convertFromVcfHeader(vcfPath.toString(),iter.getHeader());
				while(iter.hasNext())
					{
					final VariantContext ctx=iter.next();
					if(this.onlySnp && ctx.getLengthOnReference()!=1) continue;
					for(final SimplePosition m2: mutlist)
						{
						if(!m2.overlaps(ctx)) continue;
				    	if(this.onlySnp  && (ctx.getStart()!=m2.getPosition() || ctx.getEnd()!=m2.getPosition())) continue;
						report(vcfPathString,header,ctx,m2);
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
					final SimplePosition m= new SimplePosition(line);
					this.mutations.add(m);
					}
				r.close();
				}
			
			for(final String s:this.positionsList)
				{
				final SimplePosition m= new SimplePosition(s);
				this.mutations.add(m);
				}			
			
			this.out=super.openPathOrStdoutAsPrintWriter(this.outputFile);
			this.out.println("#FILE\tCHROM\tstart\tend\tID\tREF\tsample\ttype\tALLELES\tAD\tDP\tGQ");
			
			
			
			if(args.isEmpty())
				{
				try(BufferedReader br= new BufferedReader(new InputStreamReader(stdin()))) {
					br.lines().forEach(L->scanLine(L));
					}
				}
			else if(args.size()==1 && args.get(0).endsWith(".list"))
				{
				try(BufferedReader br= IOUtils.openPathForBufferedReading(Paths.get(args.get(0)))) {
					br.lines().forEach(L->scanLine(L));
					}
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
