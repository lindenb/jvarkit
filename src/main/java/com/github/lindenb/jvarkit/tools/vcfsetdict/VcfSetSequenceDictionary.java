/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfsetdict;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

The tool will try to convert the contig names ('1' -> 'chr1') according to the new dictionary.
It also fixes the contig names in the BND ALT alleles

## Example

```
java  -jar jvarkit-git/jvarkit.jar vcfsetdict --onNotFound SKIP -r ref.fasta input.vcf > out.vcf
```


END_DOC

*/
@Program(name="vcfsetdict",
	description="Set the `##contig` lines in a VCF header on the fly, fix also ",
	keywords={"vcf","dict","fai"},
	creationDate="20251014",
	modificationDate="20260225",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfSetSequenceDictionary extends OnePassVcfLauncher {
	private static final Logger LOG=Logger.of(VcfSetSequenceDictionary.class);
	private  enum OnNotFound{RAISE_EXCEPTION,SKIP,RETURN_ORIGINAL};

	@Parameter(names={"-r","-R","--reference","--dict"},description=DICTIONARY_SOURCE,required=true)
	private Path faidx=null;
	@Parameter(names={"-n","--onNotFound"},description=ContigNameConverter.OPT_ON_NT_FOUND_DESC)
	private OnNotFound onContigNotFound = OnNotFound.SKIP;			
	
	private SAMSequenceDictionary dict=null;
	
	
	public VcfSetSequenceDictionary()
		{
		}
	private String convert(final ContigNameConverter contigNameConverter, final String srcContig,final Set<String> inputContigsNotFound) {
		 final String newContig=contigNameConverter.apply(srcContig);
	
		if(StringUtils.isBlank(newContig))
			{
			if(this.onContigNotFound.equals(OnNotFound.RAISE_EXCEPTION))
				{
				throw new IllegalArgumentException("cannot convert contig "+srcContig + " for new dictionary");
				}
			else if(this.onContigNotFound.equals(OnNotFound.RETURN_ORIGINAL))
				{
				return srcContig;
				}
			else
				{
				if(!inputContigsNotFound.contains(srcContig)) {
					LOG.info("cannot convert contig "+ srcContig + " for new dictionary");
					inputContigsNotFound.add(srcContig);
					}
				return null;
				}
			}
		return newContig;
		}
	
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator in,
		final VariantContextWriter w
		) 
	    {
		final Set<String> inputContigsNotFound = new HashSet<>();
		final VCFHeader header = in.getHeader();
		final boolean has_sv_type = header.getInfoHeaderLine(VCFConstants.SVTYPE)!=null;
		final VCFHeader header2 = new VCFHeader(header);
		final ContigNameConverter contigNameConverter;
		final SAMSequenceDictionary oldDict = header.getSequenceDictionary();
		header2.setSequenceDictionary(this.dict);
		if(oldDict!=null && !oldDict.isEmpty())
			{
			contigNameConverter = ContigNameConverter.fromDictionaries(oldDict, this.dict);
			}
		else
			{
			contigNameConverter = ContigNameConverter.fromOneDictionary(this.dict);
			}
		
		w.writeHeader(header2);
		
		// for BND change:
		List<Genotype> new_genotypes  = null;
		List<Allele> new_alleles = null;

		while(in.hasNext())
			{
			final VariantContext ctx = in.next();
			boolean ok_bnd = true;
			final String newContig = convert(contigNameConverter,ctx.getContig(),inputContigsNotFound);
			// conversion failed, continue
			if(StringUtils.isBlank(newContig)) continue;
			
			//reset
			new_genotypes = null;
			new_alleles = null;
			
			if(has_sv_type && ctx.hasAttribute(VCFConstants.SVTYPE)) {
				final String svType= ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
				// if it's a breaking end, loop over the ALT alleles and try to change it
				if(svType.equals("BND")) {
					for(int i=1;i< ctx.getAlleles().size();i++)
						{
						Allele alt = ctx.getAlleles().get(i);
						if(!alt.isSymbolic()) continue;
						if(!alt.isBreakpoint()) continue;
						final String display = alt.getDisplayString();
						int colon = display.indexOf(':');
						if(colon==-1) continue;
						final char delim = display.contains("[")?'[':']';
						final int i1 = display.indexOf(delim);
						if(i1>colon) continue;
						final int i2 = colon+1==display.length()?-1:display.indexOf(delim,colon+1);
						if(!(i1+1 <colon && colon+1 <i2) ) {
							continue;
							}
						final String src_alt_contig = display.substring(i1+1,colon);
						final String bndCtg = convert(contigNameConverter,src_alt_contig,inputContigsNotFound);
						if(StringUtils.isBlank(bndCtg)) {
							ok_bnd = false;
							break;
							}
						// no conversion needed
						if(src_alt_contig.equals(bndCtg)) {
							continue;
							}
						if(ctx.hasGenotypes() && new_genotypes==null) {
							new_genotypes = new ArrayList<>(ctx.getGenotypes());
							}
						
						if(new_alleles==null) {
							new_alleles = new ArrayList<>(ctx.getAlleles());
							}
						
						
						final Allele newAlt= Allele.create(
							display.substring(0,i1+1)+
							bndCtg + 
							display.substring(colon)
							,false);
						new_alleles.set(i, newAlt);
						for(int x=0;new_genotypes!=null && x< new_genotypes.size();++x) {
							final Genotype gt = new_genotypes.get(x);
							final List<Allele> getAlleles= new ArrayList<>(gt.getAlleles());
							getAlleles.replaceAll(A->A.equals(alt)?newAlt:A);
							new_genotypes.set(x, new GenotypeBuilder(gt).alleles(getAlleles).make());
							}
						}
					}
				}
			// BND conversion failed and SKIP variant set
			if(!ok_bnd) continue;
				
			if(newContig.equals(ctx.getContig()) && new_genotypes==null && new_alleles==null) // no BND modification
				{
				w.add(ctx);
				}
			else
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.chr(newContig);
				
				if(new_alleles!=null) {
					vcb.alleles(new_alleles);
					}
				
				if(new_genotypes!=null) {
					vcb.genotypes(new_genotypes);
					}
				
				w.add(vcb.make());
				}
			}
		inputContigsNotFound.stream().forEach(chrom->
			{
			LOG.warn("Variant(s) with Contig \'"+chrom+"\' could not be converted to new Dictionary and where ignored");
			});
		inputContigsNotFound.clear();
		return 0;
	    }

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeVcf() {
		if(this.faidx==null) {
			LOG.error("no dict source defined");
			return -1;
			}
		try {
			this.dict = new SequenceDictionaryExtractor().extractRequiredDictionary(this.faidx);
			return 0;
			}
		catch (final Throwable err2)
			{
			LOG.error(err2);
			return -1;
			} 
		}

	public static void main(final String[] args) {
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
	}
}
