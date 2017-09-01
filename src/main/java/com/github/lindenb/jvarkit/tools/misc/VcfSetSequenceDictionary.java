/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
/**
BEGIN_DOC

The tool will try to convert the contig names ('1' -> 'chr1') according to the new dictionary.

## Example

```
java  -jar jvarkit-git/vcfsetdict.jar --onNotFound SKIP -r ref.fasta input.vcf > out.vcf
```


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


	@Parameter(names={"-r","-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private File faidx=null;
	@Parameter(names={"--onNotFound"},description=ContigNameConverter.OPT_ON_NT_FOUND_DESC)
	private ContigNameConverter.OnNotFound onContigNotFound =ContigNameConverter.OnNotFound.SKIP;
	
	
	@Parameter(names={"-d","--newdict"},description="At the end, save an alternate dict in that file.")
	private File newDictOut=null;

	private SAMSequenceDictionary dict=null;
	private LinkedHashMap<String, Integer> buildNewDictionary =null;

	private VcfSetSequenceDictionary()
	{
	}

	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VcfIterator in,
		final VariantContextWriter out
		) 
	    {
		final VCFHeader header=in.getHeader();		
		final VCFHeader header2=new VCFHeader(header);
		final ContigNameConverter contigNameConverter;
		
		if(this.dict!=null)
			{	
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
			}
		else
			{
			header2.setSequenceDictionary(new SAMSequenceDictionary());//
			LOG.warning("No sequence dictionary was defined");
			contigNameConverter = ContigNameConverter.getIdentity();
			}
		
		contigNameConverter.setOnNotFound(this.onContigNotFound);
		
		final Set<String> inputContigsNotFound=new HashSet<>();
		
		out.writeHeader(header2);
		final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(header).logger(LOG);
		while(in.hasNext())
		{
			final VariantContext ctx=progress.watch(in.next());
			
			
			if(this.buildNewDictionary!=null)
			{
				Integer length=this.buildNewDictionary.get(ctx.getContig());
				if(length==null) length=0;
				if(ctx.getEnd()>length)
				{
					this.buildNewDictionary.put(ctx.getContig(),ctx.getEnd());
				}
			}

			final String newContig = contigNameConverter.apply(ctx.getContig());
			if(newContig==null)
				{
				if(!inputContigsNotFound.contains(ctx.getContig())) {
					LOG.info("cannot convert contig "+ctx.getContig()+ " for new dictionary");
					inputContigsNotFound.add(ctx.getContig());
					}
				continue;
				}
			else if(newContig.equals(ctx.getContig()))
				{
				out.add(ctx);
				}
			else
				{
				out.add(new VariantContextBuilder(ctx).chr(newContig).make());
				}
			}
		progress.finish();
		
		if(!inputContigsNotFound.isEmpty())
			{
			for(final String chrom: inputContigsNotFound)
				{
				LOG.warn("Variant with Contig "+chrom+" could not be converted to new Dictionary and where ignored");
				}
			}
		
		return 0;
	}

	@Override
	public int doWork(final List<String> args) {
		if (newDictOut != null) {
			if (!newDictOut.getName().endsWith(".dict")) {
				LOG.error("dictionary should end with .dict :" + newDictOut);
				return -1;
			}
			this.buildNewDictionary = new LinkedHashMap<String, Integer>();
		}
		FileWriter out = null;
		try {
			if (this.faidx != null) {
				this.dict = SAMSequenceDictionaryExtractor.extractDictionary(faidx);
			} else 
				if(newDictOut==null)
				{
				LOG.error("new dictionary undefined");
				return -1;
				}

			final int err = doVcfToVcf(args, this.outputFile);

			if (newDictOut != null) {
				LOG.info("Saving alt dictionary " + newDictOut);

				final List<SAMSequenceRecord> list = new ArrayList<SAMSequenceRecord>(buildNewDictionary.size());
				for (final String k : this.buildNewDictionary.keySet()) {
					list.add(new SAMSequenceRecord(k, this.buildNewDictionary.get(k)));
				}
				final SAMFileHeader sfh = new SAMFileHeader();
				sfh.setSequenceDictionary(new SAMSequenceDictionary(list));
				final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
				codec.setValidationStringency(htsjdk.samtools.ValidationStringency.SILENT);
				out = new FileWriter(this.newDictOut);
				codec.encode(out, sfh);
				out.flush();
			}
			
			return err;
		} catch (final Exception err2) {
			LOG.error(err2);
			return -1;
		} finally {
			CloserUtil.close(out);
		}

	}

	public static void main(final String[] args)
	{
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
	}
}
