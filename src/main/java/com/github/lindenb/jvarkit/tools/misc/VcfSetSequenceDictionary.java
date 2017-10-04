/*
The MIT License (MIT)

Copyright (c) 2017 Pierre Lindenbaum

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

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlType;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VariantContextWriterFactory;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
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
	keywords={"vcf","dict","fai"}
	)
public class VcfSetSequenceDictionary extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfSetSequenceDictionary.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@ParametersDelegate
	private CtxWriterFactory component = new CtxWriterFactory();
	

	@XmlType(name="vcfsetdict")
	@XmlAccessorType(XmlAccessType.FIELD)
	public static class CtxWriterFactory 
		implements VariantContextWriterFactory
			{
			@XmlElement(name="reference")
			@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
			private File faidx=null;
			
			@XmlElement(name="on-not-found")
			@Parameter(names={"--onNotFound"},description=ContigNameConverter.OPT_ON_NT_FOUND_DESC)
			private ContigNameConverter.OnNotFound onContigNotFound =ContigNameConverter.OnNotFound.SKIP;			
		
			@XmlTransient
			private SAMSequenceDictionary dict=null;
			
			public void setReference(final File faidx) {
				this.faidx = faidx;
				}
			public void setOnContigNotFound(final ContigNameConverter.OnNotFound onContigNotFound) {
				this.onContigNotFound = onContigNotFound;
				}
			
			
			
			private  class CtxWriter extends DelegateVariantContextWriter
				{
				private ContigNameConverter contigNameConverter;
				private final Set<String> inputContigsNotFound=new HashSet<>();

				CtxWriter(final VariantContextWriter delegate) {
					super(delegate);
					}
				
				@Override
				public void writeHeader(final VCFHeader header) {
					final VCFHeader header2 = new VCFHeader(header);
					
					final SAMSequenceDictionary oldDict = header.getSequenceDictionary();
					header2.setSequenceDictionary(CtxWriterFactory.this.dict);
					if(oldDict!=null && !oldDict.isEmpty())
						{
						this.contigNameConverter = ContigNameConverter.fromDictionaries(oldDict, CtxWriterFactory.this.dict);
						}
					else
						{
						this.contigNameConverter = ContigNameConverter.fromOneDictionary(CtxWriterFactory.this.dict);
						}
					this.contigNameConverter.setOnNotFound(CtxWriterFactory.this.onContigNotFound);
					super.writeHeader(header2);
					}
				@Override
				public void add(final VariantContext ctx) {
					final String newContig = this.contigNameConverter.apply(ctx.getContig());
					if(newContig==null)
						{
						if(!inputContigsNotFound.contains(ctx.getContig())) {
							LOG.info("cannot convert contig "+ctx.getContig()+ " for new dictionary");
							inputContigsNotFound.add(ctx.getContig());
							}
						return;
						}
					else if(newContig.equals(ctx.getContig()))
						{
						super.add(ctx);
						}
					else
						{
						super.add(new VariantContextBuilder(ctx).chr(newContig).make());
						}
					}
				@Override
				public void close() {
					this.inputContigsNotFound.stream().forEach(chrom->
						{
						LOG.warn("Variant(s) with Contig \'"+chrom+"\' could not be converted to new Dictionary and where ignored");
						});
					this.inputContigsNotFound.clear();
					super.close();
					}
				}
			
			@Override
			public int initialize() {
				Objects.requireNonNull(this.faidx);
				this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidx);
				return 0;
				}
			
			@Override
			public VariantContextWriter open(VariantContextWriter delegate) {
				return new CtxWriter(delegate);
				}
			
			@Override
			public void close() throws IOException {
				this.dict=null;
				}
			
			}

	public VcfSetSequenceDictionary()
		{
		}

	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VcfIterator in,
		final VariantContextWriter delegate
		) 
	    {
		final VariantContextWriter out = this.component.open(delegate);
		
		out.writeHeader(in.getHeader());
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(in.getHeader()).logger(LOG);
		while(in.hasNext())
			{
			out.add(progress.watch(in.next()));
			}
		progress.finish();
		return 0;
	    }

	@Override
	public int doWork(final List<String> args) {
		try {
			if(this.component.initialize()!=0) {
				return -1;
				}
			 return doVcfToVcf(args, this.outputFile);
		} catch (final Exception err2) {
			LOG.error(err2);
			return -1;
		} finally {
			CloserUtil.close(this.component);
		}
	}

	public static void main(final String[] args) {
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
	}
}
