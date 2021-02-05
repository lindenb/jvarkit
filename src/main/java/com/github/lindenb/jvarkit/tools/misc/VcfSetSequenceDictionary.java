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

import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

The tool will try to convert the contig names ('1' -> 'chr1') according to the new dictionary.

## Example

```
java  -jar jvarkit-git/vcfsetdict.jar --onNotFound SKIP -r ref.fasta input.vcf > out.vcf
```

## History

* [20170906] remove the creation of a dictionary, moved to VcfCreateDictionary


END_DOC

*/
@Program(name="vcfsetdict",
	description="Set the `##contig` lines in a VCF header on the fly",
	keywords={"vcf","dict","fai"},
	creationDate="20140105",
	modificationDate="20210201"
	)
public class VcfSetSequenceDictionary extends OnePassVcfLauncher {
	private static final Logger LOG=Logger.build(VcfSetSequenceDictionary.class).make();
	private  enum OnNotFound{RAISE_EXCEPTION,SKIP,RETURN_ORIGINAL};

	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path faidx=null;
	@Parameter(names={"-n","--onNotFound"},description=ContigNameConverter.OPT_ON_NT_FOUND_DESC)
	private OnNotFound onContigNotFound = OnNotFound.SKIP;			
	
	private SAMSequenceDictionary dict=null;
	
	
	public VcfSetSequenceDictionary()
		{
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
		
		
		while(in.hasNext())
			{
			final VariantContext ctx = in.next();
			
			String newContig = contigNameConverter.apply(ctx.getContig());
			
			
			if(StringUtils.isBlank(newContig))
				{
				if(this.onContigNotFound.equals(OnNotFound.RAISE_EXCEPTION))
					{
					LOG.error("cannot convert contig "+ctx.getContig()+ " for new dictionary");
					return -1;
					}
				else if(this.onContigNotFound.equals(OnNotFound.RETURN_ORIGINAL))
					{
					newContig = ctx.getContig();
					}
				else
					{
					if(!inputContigsNotFound.contains(ctx.getContig())) {
						LOG.info("cannot convert contig "+ctx.getContig()+ " for new dictionary");
						inputContigsNotFound.add(ctx.getContig());
						}
					continue;
					}
				}
			else if(newContig.equals(ctx.getContig()))
				{
				w.add(ctx);
				}
			else
				{
				w.add(new VariantContextBuilder(ctx).chr(newContig).make());
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
			LOG.error("REF not defined");
			return -1;
			}
		try {
			this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidx);
			return 0;
			}
		catch (final Exception err2)
			{
			LOG.error(err2);
			return -1;
			} 
		}

	public static void main(final String[] args) {
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
	}
}
