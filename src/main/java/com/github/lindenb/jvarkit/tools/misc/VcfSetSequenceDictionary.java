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
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
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
	keywords={"vcf","dict","fai"}
	)
public class VcfSetSequenceDictionary extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfSetSequenceDictionary.class).make();
	private  enum OnNotFound{RAISE_EXCEPTION,SKIP,RETURN_ORIGINAL};

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private File faidx=null;
	@Parameter(names={"-n","--onNotFound"},description=ContigNameConverter.OPT_ON_NT_FOUND_DESC)
	private OnNotFound onContigNotFound = OnNotFound.SKIP;			
	@Parameter(names={"-ho","--header-only"},description=
			"only change the vcf header. Keep the whole VCF body unchanged. The idea is to use a faster(?) `sed sed 's/^chr//' ` for the VCF body. " )
	private boolean header_only=false;			

	
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
		
		
		final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().
				dictionary(header).
				validatingContig(false).
				logger(LOG).
				build();
		while(in.hasNext())
			{
			final VariantContext ctx=progress.apply(in.next());
			
			String newContig = contigNameConverter.apply(ctx.getContig());
			
			
			if(newContig==null)
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
		progress.close();
		inputContigsNotFound.stream().forEach(chrom->
			{
			LOG.warn("Variant(s) with Contig \'"+chrom+"\' could not be converted to new Dictionary and where ignored");
			});
		inputContigsNotFound.clear();
		return 0;
	    }

	private int headerOnly(final String inputName,final File outputFile) throws IOException {
		BufferedReader br = null;
		PrintWriter pw = null;
		try {
			br = super.openBufferedReader(inputName);
			pw = super.openFileOrStdoutAsPrintWriter(outputFile);
			VCFUtils.CodecAndHeader cah =  VCFUtils.parseHeader(br);
			final VCFHeader header2 = new VCFHeader(cah.header);
			header2.setSequenceDictionary(this.dict);
			final ByteArrayOutputStream baos=new ByteArrayOutputStream();
			final VariantContextWriter vcw=VCFUtils.createVariantContextWriterToOutputStream(baos);
			vcw.writeHeader(header2);
			vcw.close();
			pw.print(new String(baos.toByteArray()));
			IOUtils.copyTo(br, pw);
			pw.flush();
			pw.close();pw=null;
			br.close();br=null;
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
		} finally {
			CloserUtil.close(br);
			CloserUtil.close(pw);
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.faidx==null) {
			LOG.error("REF not defined");
			return -1;
			}
		try {
			this.dict = SAMSequenceDictionaryExtractor.extractDictionary(this.faidx);
			if(this.header_only)
				{
				return headerOnly(oneFileOrNull(args), this.outputFile);
				}
			else
				{
				return doVcfToVcf(args, this.outputFile);
				}
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
