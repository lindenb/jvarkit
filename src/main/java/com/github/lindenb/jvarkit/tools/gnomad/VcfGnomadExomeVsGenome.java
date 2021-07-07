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
package com.github.lindenb.jvarkit.tools.gnomad;

import java.nio.file.Path;
import java.util.Collections;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.vcf.BufferedVCFReader;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;

@Program(name="vcfgnomadexomevsgenome",
description="filter out variant in exome but not in genome",
keywords={"vcf","annotation","gnomad"},
modificationDate="20210707",
creationDate="20210707"
)
public class VcfGnomadExomeVsGenome extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.build(VcfGnomadExomeVsGenome.class).make();
	@Parameter(names={"-exome"},description="Path to Indexed Gnomad VCF exome file.",required=true)
	private Path exomePath =null;
	@Parameter(names={"-genome"},description="Path to Indexed Gnomad VCF genome file.",required=true)
	private Path genomePath =null;
	@Parameter(names={"--bufferSize"},description= BufferedVCFReader.OPT_BUFFER_DESC+" "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int gnomadBufferSize= 10_000;
	@Parameter(names={"--filter"},description= BufferedVCFReader.OPT_BUFFER_DESC+" "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private String filterSr= "GNOMAD_ONLY_EXOME";

	private BufferedVCFReader exomeReader = null;
	private BufferedVCFReader genomeReader = null;
	private ContigNameConverter ctgNameConverter=null; 
	
	@Override
	protected int beforeVcf()
		{
		final UnaryOperator<VariantContext> simplifier = V-> new VariantContextBuilder(V).noGenotypes().noID().attributes(Collections.emptyMap()).make();
		SAMSequenceDictionary dict1=null;
		try {
			final VCFReader r = VCFReaderFactory.makeDefault().open(this.exomePath,true);
			this.exomeReader = new BufferedVCFReader(r, this.gnomadBufferSize);
			this.exomeReader.setSimplifier(simplifier);
			this.ctgNameConverter = ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(r.getHeader()));
			dict1 = SequenceDictionaryUtils.extractRequired(r.getHeader());
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		

		try {
			final VCFReader r = VCFReaderFactory.makeDefault().open(this.genomePath,true);
			this.genomeReader = new BufferedVCFReader(r, this.gnomadBufferSize);
			this.genomeReader.setSimplifier(simplifier);
			final SAMSequenceDictionary dict2 = SequenceDictionaryUtils.extractRequired(r.getHeader());
			SequenceUtil.assertSequenceDictionariesEqual(dict1, dict2);
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		return super.beforeVcf();
		}
	
	boolean isVariantIn(final VariantContext ctx,final BufferedVCFReader reader) {
			final String ctg = this.ctgNameConverter.apply(ctx.getContig());
			if(StringUtils.isBlank(ctg)) return false;
			try(CloseableIterator<VariantContext> iterE = reader.query(new SimpleInterval(ctg,ctx.getStart(),ctx.getEnd()))) {
				while(iterE.hasNext()) {
					final VariantContext ctx2 = iterE.next();
					if(ctx.getStart()!=ctx2.getStart()) continue;
					if(!ctx.getReference().equals(ctx2.getReference())) continue;
					final Set<Allele> alt2 = ctx2.getAlternateAlleles().stream().filter(A->AcidNucleics.isATGC(A)).collect(Collectors.toSet());

					if(ctx.getAlternateAlleles().stream().anyMatch(A->alt2.contains(A))) return true;
					}
				}
			return false;
			}

	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin,
			VariantContextWriter out)
		{
		final VCFHeader h2 = new VCFHeader(iterin.getHeader());
		long removed = 0L;
		if(!StringUtils.isBlank(this.filterSr)) {
			final VCFFilterHeaderLine fh2 = new VCFFilterHeaderLine(
					this.filterSr,
					"variant wasfound in "+this.exomePath+" but not in "+this.genomePath
					);
			h2.addMetaDataLine(fh2);
			}
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		while(iterin.hasNext()) {
			final VariantContext ctx = iterin.next();
			if(isVariantIn(ctx,this.exomeReader) && !isVariantIn(ctx,this.genomeReader)) {
				++removed;
				if(StringUtils.isBlank(this.filterSr)) {
					continue;
					}
				else
					{
					out.add(new VariantContextBuilder(ctx).filter(this.filterSr).make());
					}
				}
			else
				{
				out.add(ctx);
				}
			}
		if(removed>0L) LOG.info(String.valueOf(removed)+" variants were found in "+this.exomePath+" but not in "+this.genomePath);
		return 0;
		}
	
	@Override
	protected void afterVcf() {
		try {
			this.exomeReader.close();
			}
		catch(final Throwable err) {
			LOG.error(err);
			}
		try {
			this.genomeReader.close();
			}
		catch(final Throwable err) {
			LOG.error(err);
			}
		}
	
	public static void main(final String[] args)
		{
		new VcfGnomadExomeVsGenome().instanceMainWithExit(args);

		}

	}
