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
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;


import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

@Program(name="vcfsetdict",description="Set the ##contig lines in a VCF header",
	deprecatedMsg="Use picard UpdateVcfSequenceDictionary",
	keywords={"vcf","dict","fai"}
	)
public class VcfSetSequenceDictionary extends Launcher
{
	private static final Logger LOG=Logger.build(VcfSetSequenceDictionary.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;


	@Parameter(names={"-r","-R","--reference"},description="indexed reference")
	private File faidx=null;

	@Parameter(names={"-d"},description="at the end, save an alternate dict in that file.")
	private File newDictOut=null;

	private SAMSequenceDictionary dict=null;
	private LinkedHashMap<String, Integer> newdict=null;

	private VcfSetSequenceDictionary()
	{
	}

	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final VCFHeader header=in.getHeader();
		final Set<VCFHeaderLine> meta2=new LinkedHashSet<VCFHeaderLine>();
		for(final VCFHeaderLine L:header.getMetaDataInInputOrder())
		{
			if(!L.getKey().equals(VCFHeader.CONTIG_KEY))
			{
				meta2.add(L);
			}
		}

		meta2.add(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		meta2.add(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));

		if(dict!=null)
		{	
			meta2.addAll(VCFUtils.samSequenceDictToVCFContigHeaderLine(dict));
		}
		else
		{
			LOG.warning("No sequence dictionary was defined");
		}
		final VCFHeader header2=new VCFHeader(meta2, header.getSampleNamesInOrder());
		out.writeHeader(header2);

		while(in.hasNext())
		{
			final VariantContext ctx=in.next();
			if(dict!=null && dict.getSequenceIndex(ctx.getContig())==-1)
			{
				LOG.warning("Unknown chromosome "+ctx.getContig());
			}
			if(newdict!=null)
			{
				Integer length=this.newdict.get(ctx.getContig());
				if(length==null) length=0;
				if(ctx.getEnd()>length)
				{
					this.newdict.put(ctx.getContig(),ctx.getEnd());
				}
			}


			out.add(ctx);
		}
		return 0;
	}

	@Override
	public int doWork(List<String> args) {
		if (newDictOut != null) {
			if (!newDictOut.getName().endsWith(".dict")) {
				LOG.error("dictionary should end with .dict :" + newDictOut);
				return -1;
			}
			this.newdict = new LinkedHashMap<String, Integer>();
		}
		FileWriter out = null;
		try {
			if (this.faidx != null) {
				this.dict = SAMSequenceDictionaryExtractor.extractDictionary(faidx);
			}

			int err = doVcfToVcf(args, outputFile);

			if (newDictOut != null) {
				LOG.info("Saving alt dictionary " + newDictOut);

				final List<SAMSequenceRecord> list = new ArrayList<SAMSequenceRecord>(newdict.size());
				for (String k : this.newdict.keySet()) {
					list.add(new SAMSequenceRecord(k, this.newdict.get(k)));
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
		} catch (Exception err2) {
			LOG.error(err2);
			return -1;
		} finally {
			CloserUtil.close(out);
		}

	}

	public static void main(String[] args)
	{
		new VcfSetSequenceDictionary().instanceMainWithExit(args);
	}
}
