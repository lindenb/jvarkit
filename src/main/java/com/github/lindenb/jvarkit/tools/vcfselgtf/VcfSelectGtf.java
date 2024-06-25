/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfselgtf;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.MinMaxInteger;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.tools.jvarkit.JvarkitCentral;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFCodec;
import com.github.lindenb.jvarkit.util.bio.gtf.GTFLine;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC


## Motivation

select variants in a vcf according to their position in a transcript


variant must be annotated with snpeff or vep

END_DOC
*/
@Program(name="vcfselectgtf",
	description="Select variant by their position in the transcripts in a GTF (intron, first-intron, utr, etc...)",
	keywords={"vcf","gtf"},
	creationDate="20240725",
	modificationDate="20240725",
	jvarkit_amalgamion =  true,
	menu="VCF Manipulation"
	)
public class VcfSelectGtf extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.build(VcfSelectGtf.class).make();
	private enum InfoType {GENE_ID,TRANSCRIPT_ID};
	private enum Coding {all,protein_coding};
	@Parameter(names={"-gtf","--gtf"},description="GTF file indexed with tabix",required = true)
	private String gtfPath;

	private TabixReader tabixReader=null;
	private Interval bufferInterval = null;
	private List<GTFLine> buffer=new ArrayList<>();
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	private class GTFNode
		{
		final GTFLine line;
		final List<GTFNode> children=new ArrayList<>();
		GTFNode(final GTFLine line) {
			this.line = line;
			}
		}
	
	@Override
	protected int beforeVcf() {
		try {
			this.tabixReader = new TabixReader(this.gtfPath);
			}
		catch (final IOException err) {
			LOG.error(err);
			return -1;
			}
		
		return super.beforeVcf();
		}
	
	@Override
	protected void afterVcf() {
		if(this.tabixReader!=null) this.tabixReader.close();
		this.tabixReader=null;
		super.afterVcf();
		}
	
	private void testVariant(VariantContext ctx, Set<String> transctriptIds) throws IOException {
		final GTFCodec gtfCodec = new GTFCodec();

		final MinMaxInteger minMaxPos = new MinMaxInteger();
		final Map<String,GTFNode> genes = new HashMap<>();
		TabixReader.Iterator iter= tabixReader.query(ctx.getContig()+":"+ctx.getStart()+"-"+ctx.getEnd());
		for(;;) {
			String line = iter.next();
			if(line==null) break;
			final GTFLine  gtf= gtfCodec.decode(line);
			if(!gtf.getType().equals("gene")) continue;
			if(StringUtils.isBlank(gtf.getGeneId())) continue;
			minMaxPos.accept(gtf.getStart());
			minMaxPos.accept(gtf.getEnd());
			GTFNode node= new GTFNode(gtf);
			genes.put(gtf.getGeneId(), node);
			}
		if(!minMaxPos.isEmpty()) {
			iter= tabixReader.query(ctx.getContig()+":"+minMaxPos.getMinAsInt()+"-"+minMaxPos.getMaxAsInt());
			for(;;) {
				String line = iter.next();
				if(line==null) break;
				final GTFLine  gtf= gtfCodec.decode(line);
				if(gtf.getType().equals("gene")) continue;
				if(StringUtils.isBlank(gtf.getGeneId())) continue;
				GTFNode node=genes.get(gtf.getGeneId());
				if(node==null) continue;
				}
			}
		}
	
	@Override
	protected int doVcfToVcf(final String inputName, VCFIterator iterin, VariantContextWriter out) {
		try {
			final VCFHeader header = iterin.getHeader();
			
			final VepPredictionParser vep = new VepPredictionParserFactory(header).get();
			final AnnPredictionParser ann = new AnnPredictionParserFactory(header).get();
			if(!vep.isValid() && !ann.isValid()) {
				LOG.error("no valid functional prediction found. Did you annotate the VCF with vep and/or snpeff ?");
				return -1;
				}
			final Function<VariantContext,Set<String>> vc2transcript = CTX->{
				final Set<String> set = new HashSet<>();
				if(vep.isValid()) {
					vep.getPredictions(CTX).
						stream().
						map(PRED->PRED.getFeature()).
						filter(S->!StringUtils.isBlank(S)).
						forEach(S->set.add(S));
					}
				if(ann.isValid()) {
					ann.getPredictions(CTX).
						stream().
						map(PRED->PRED.getFeatureId()).
						filter(S->!StringUtils.isBlank(S)).
						forEach(S->set.add(S));
					}
				return set;
				};
			
			JVarkitVersion.getInstance().addMetaData(this, header);
			
			out.writeHeader(header);
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				final Set<String> transcripts_id= vc2transcript.apply(ctx);
				
				
				out.add(ctx);
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(final String[] args)
		{
		new VcfSelectGtf().instanceMainWithExit(args);
		}

}
