/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.findavariation;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFReader;
/** 
BEGIN_DOC
 
## Example

```
$ find ./ -name "*.vcf" -o -name "*.vcf.gz" |\
   java -jar dist/findavariation.jar -V search.vcf


htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12878	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12891	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	1	8216713	8216713	yossi-1		NA12892	HET	A G
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12878	HOM_REF	C C
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12891	HET	C T
htsjdk/testdata/htsjdk/samtools/intervallist/IntervalListFromVCFTestManual.vcf	2	2	2	.		NA12892	HET	C T
```

 
END_DOC
*/
@Program(name="findavariation",
	description="Finds a specific mutation in a list of VCF files",
	keywords={"vcf","variation","search","find","bcf"},
	creationDate="20140623",
	modificationDate="202050902",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class FindAVariation extends Launcher
	{
	private enum Matcher {overlap,id,chrom_pos,chrom_pos_ref,chrom_pos_any_alt,chrom_pos_all_alt};
	private static final Logger LOG = Logger.of(FindAVariation.class);
	@Parameter(names={"-V","--vcf","--variants"}, description="A VCF file containing the variants",required = true)
	private Path vcfFile = null;
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-homref","--homref"},description="Hide HOM_REF genotypes")
	private boolean hideHomRef=false;
	@Parameter(names={"-nocall","--nocall"},description="Hide NO_CALL genotypes")
	private boolean hideNoCall=false;
	@Parameter(names={"--enable-no-index"},description="Enabled scan of non-indexed files.")
	private boolean enabled_non_indexed=false;
	@Parameter(names={"--pair-logic","-p"},description=" Matching records by 'x'")
	private Matcher matcher= Matcher.chrom_pos_any_alt;

	
	
    public FindAVariation()
    	{
    	}		
   
    
    private static boolean overlap(final Matcher matcher,final VariantContext userVar, final VariantContext vcfVar) {
    	if(!userVar.overlaps(vcfVar)) return false;
    	switch(matcher) {
    		case overlap: return true;
    		case id:
    			{
				if(!userVar.hasID() || !vcfVar.hasID() ) return false;
				final Set<String> set = Arrays.stream( CharSplitter.COMMA.split(vcfVar.getID())).collect(Collectors.toSet());
				for(String id : CharSplitter.COMMA.split(userVar.getID())) {
					if(id.equals(VCFConstants.EMPTY_ID_FIELD)) continue;
					if(set.contains(id)) return true;
					}
				return false;
    			}
    		case chrom_pos:
    			return userVar.contigsMatch(vcfVar) && userVar.getStart() == vcfVar.getStart();
    		case chrom_pos_ref:
    			return overlap(Matcher.chrom_pos,userVar,vcfVar) && userVar.getNAlleles()>0 && userVar.getReference().equals(vcfVar.getReference());
    		case chrom_pos_any_alt:
    			{
    			if(!overlap(Matcher.chrom_pos_ref,userVar,vcfVar)) return false;
    			final Set<Allele> set =new HashSet<>(userVar.getAlternateAlleles());
    			for(Allele a: vcfVar.getAlternateAlleles()) {
    				if(set.contains(a)) return true;
    				}
    			return false;
    			}
    		case chrom_pos_all_alt:
				{
				if(!overlap(Matcher.chrom_pos_ref,userVar,vcfVar)) return false;
				final Set<Allele> set1 =new HashSet<>(userVar.getAlternateAlleles());
				final Set<Allele> set2 =new HashSet<>(vcfVar.getAlternateAlleles());
				return set1.containsAll(set2);
				}
   			}
    	return false;
    	}
    
    private boolean overlap(final VariantContext userVar, final VariantContext vcfVar) {
    	return overlap(this.matcher,userVar,vcfVar);
    	}
    
    private void reportPos(	final PrintWriter out,final String f,final VCFHeader header,final VariantContext ctx)
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
    		final PrintWriter out,
    		final String f,
    		final VCFHeader header,
    		final VariantContext ctx,
    		final VariantContext mpos
    		)
    	{
    	if(!ctx.hasGenotypes())
    		{
    		reportPos(out,f,header,ctx);
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
    			reportPos(out,f,header,ctx);
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
    
    private List<VariantContext> convertFromVcfHeader(final String f,final VCFHeader h,final List<VariantContext> mutations)
    	{
    	final SAMSequenceDictionary dict= h.getSequenceDictionary();
    	if(dict==null || dict.isEmpty()) {
    		return mutations;
    		}
    	
    	final List<VariantContext> copy=new ArrayList<>(mutations.size());
    	final ContigNameConverter ctgCvt = ContigNameConverter.fromOneDictionary(dict);
    	
    	for(final VariantContext m:mutations)
    		{
    		final String s=ctgCvt.apply(m.getContig());
    		if(StringUtils.isBlank(s))
    			{
    			LOG.warn("Cannot convert chrom "+m.getContig()+" in "+f);
    			continue;
    			}
    		copy.add(new VariantContextBuilder(m).chr(s).make());
    		}
    	return copy;
    	}

    
    
    private void scanLine(final PrintWriter out,final String line,final List<VariantContext> mutations,VCFHeader userHeader)  {
    	if(StringUtils.isBlank(line) || line.startsWith("#")) return;
    	if(IOUtil.isUrl(line)) {
    		scanRemote(out,line,mutations,userHeader);
    		}
    	else
    		{
    		scanPath(out,line,mutations,userHeader);
    		}
    	}
    
    
    private void scanRemote(final PrintWriter out,final String url,final List<VariantContext> mutations,final VCFHeader userHeader)  {
    	try (VCFReader tabix=VCFReaderFactory.makeDefault().open(url)) {
    		final VCFHeader header = tabix.getHeader();
    		final SAMSequenceDictionary dict = header.getSequenceDictionary();
    		for(final VariantContext m:convertFromVcfHeader(url,header,mutations))
    			{
    			if(dict.getSequence(m.getContig())==null) continue;
    			final java.util.Iterator<VariantContext> iter2 = tabix.query(m);
    		    while(iter2.hasNext())
                      {
                      final VariantContext ctx=iter2.next();
                      if(!overlap(m,ctx)) continue;
                      report(out,url,header,ctx,m);
                      }
                  }
    		}
		catch(final Throwable err)
			{
			LOG.severe("URL= "+url,err);
			}
    	}
    
    private void scanPath(final PrintWriter out,final String vcfPathString,final List<VariantContext> mutations)  {
    	final Path vcfPath=Paths.get(vcfPathString);
    	if(Files.isDirectory(vcfPath)) return;
    	if(!Files.isReadable(vcfPath)) return;
    	Throwable lastError=null;

		for(int n_try=0;n_try<2;++n_try) {
			lastError=null;
			try(VCFReader r= VCFReaderFactory.makeDefault().open(vcfPath,n_try==0)) {
				final VCFHeader header =r.getHeader();
				if(n_try==0) {
					for(final VariantContext m : convertFromVcfHeader(vcfPath.toString(),header, mutations))
						{
						try( CloseableIterator<VariantContext> iter2 = r.query(m)) {
							while(iter2.hasNext())
								{
								final VariantContext ctx=iter2.next();
			                    if(!overlap(m,ctx)) continue;
								report(out,vcfPathString,header,ctx,m);
								}
							}
						}
					}
				else
					{
					final List<VariantContext> mutlist = convertFromVcfHeader(vcfPath.toString(),r.getHeader(), mutations);
					try( CloseableIterator<VariantContext> iter = r.iterator()) {
						while(iter.hasNext())
							{
							final VariantContext ctx=iter.next();
							for(final VariantContext m: mutlist)
								{
								if(!overlap(m,ctx)) continue;
								report(out,vcfPathString,header,ctx,m);
								}
							}
						}
					}
				}   	
			catch(Throwable err) {
				lastError=err;
				}
			
			if(!this.enabled_non_indexed) break;
			}
		if(lastError!=null) {
			LOG.warn(vcfPathString+" "+lastError.getMessage());
			}
		}
    

	
	@Override
	public int doWork(final List<String> args) {
		try
			{
			final List<VariantContext> mutations=new ArrayList<>();
			VCFHeader userHeader=null;
			try(final VCFReader r = VCFReaderFactory.makeDefault().open(this.vcfFile,false) ) {
				userHeader = r.getHeader();
				try(CloseableIterator<VariantContext> iter=r.iterator()) {
					while(iter.hasNext()) {
						mutations.add(new VariantContextBuilder(iter.next()).noGenotypes().make());
						}
					}
				}
			
			
			try(PrintWriter out=super.openPathOrStdoutAsPrintWriter(this.outputFile) ) {
				out.println("#FILE\tCHROM\tstart\tend\tID\tREF\tsample\ttype\tALLELES\tAD\tDP\tGQ");
				
				if(args.isEmpty())
					{
					try(BufferedReader br= new BufferedReader(new InputStreamReader(stdin()))) {
						br.lines().forEach(L->scanLine(out,L,mutations,userHeader));
						}
					}
				else if(args.size()==1 && args.get(0).endsWith(".list"))
					{
					try(BufferedReader br= IOUtils.openPathForBufferedReading(Paths.get(args.get(0)))) {
						br.lines().forEach(L->scanLine(out,L,mutations,userHeader));
						}
					}
				else 
					{
					args.forEach(L->scanLine(out,L,mutations,userHeader));
					}
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
	public static void main(final String[] args) {
		new FindAVariation().instanceMainWithExit(args);

	}

}
