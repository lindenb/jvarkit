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
package com.github.lindenb.jvarkit.tools.retrocopy;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.ucsc.UcscTranscript;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptCodec;
import com.github.lindenb.jvarkit.ucsc.UcscTranscriptReader;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example

```
$ java -jar dist/knownretrocopy.jar --gtf Homo_sapiens.GRCh37.87.gtf.gz candidateSV.vcf.gz | grep RETR

##FILTER=<ID=RETROCOPY_INTRON,Description="variant could be a deleted intron from a retrocopy">
##FILTER=<ID=RETROCOPY_KNOWN,Description="variant could be a deleted intron from a known retrocopy">
##INFO=<ID=RETROCOPY,Number=.,Type=String,Description="Identifiers for the retrocopies.">
chr1	38077349	MantaDEL:2204:0:0:0:0:0	CTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTA	C	.RETROCOPY_INTRON	CIGAR=1M71D;DOWNSTREAM_PAIR_COUNT=0;END=38077420;PAIR_COUNT=0;RETROCOPY=ENSG00000169218,ENST00000356545,ENST00000401071,RSPO1,RSPO1-201,RSPO1-202;SVLEN=-71;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=0
chr1	38077349	MantaDEL:2204:0:1:0:0:0	CTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTA	C	.RETROCOPY_INTRON	CIGAR=1M71D;DOWNSTREAM_PAIR_COUNT=0;END=38077420;PAIR_COUNT=0;RETROCOPY=ENSG00000169218,ENST00000356545,ENST00000401071,RSPO1,RSPO1-201,RSPO1-202;SVLEN=-71;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=0
chr1	38077349	MantaDEL:2204:1:1:0:0:0	CTTATCCCCATACTAGTTATTATCGAAACCATCAGCCTACTCATTCAACCAATAGCCCTGGCCGTACGCCTA	C	.RETROCOPY_INTRON	CIGAR=1M71D;DOWNSTREAM_PAIR_COUNT=0;END=38077420;PAIR_COUNT=0;RETROCOPY=ENSG00000169218,ENST00000356545,ENST00000401071,RSPO1,RSPO1-201,RSPO1-202;SVLEN=-71;SVTYPE=DEL;UPSTREAM_PAIR_COUNT=0
(...)
```

END_DOC

*/
@Program(name="knownretrocopy",
description="Annotate VCF structural variants that could be intron from retrocopies.",
keywords={"gtf","retrocopy","deletion"},
creationDate="20190815",
modificationDate="20230817",
jvarkit_amalgamion = true
)
public class KnownRetroCopy extends OnePassVcfLauncher
	{
	private static final Logger LOG = Logger.of(KnownRetroCopy.class);
	@Parameter(names={"-genpred","--genpred","--kg","--transcripts"},description= UcscTranscriptReader.OPT_DESC,required=true)
	private String genPredURI = null;
	@Parameter(names={"-d","--distance"},description="max distance between an intron and the deletion found in the VCF")
	private int distance = 10;
	@Parameter(names={"--mic","--min-intron-count"},description="Min intron count.",hidden=true)
	private int min_intron_count = 0;
	@Parameter(names={"-k","--known"},description="Gene-ID of known retrogenes. One per line. A source could be : http://retrogenedb.amu.edu.pl/static/download/")
	private Path knownPath = null;
	
	private final static String ATT_RETROCOPY="RETROCOPY";
	private final static String ATT_FILTER_INTRON="RETROCOPY_INTRON";
	private final static String ATT_FILTER_KNOWN="RETROCOPY_KNOWN";

	
	private static class Intron extends SimpleInterval {
		final String transcriptId;
		final String geneId;
		Intron(UcscTranscript.Intron intron) {
			super(intron);
			this.transcriptId = intron.getTranscript().getTranscriptId();
			this.geneId = intron.getTranscript().getName2();
		}
	}
	
	private boolean isWithinDistance(int a, int b) {
		return Math.abs(a-b) < this.distance;
	}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {	
		try {
			/* load the reference genome */
			/* create a contig name converter from the REF */
	
			final Set<String> knownGeneIds;
			if(this.knownPath!=null) {
				try(BufferedReader br = IOUtils.openPathForBufferedReading(this.knownPath)) {
					knownGeneIds = br.lines().
						filter(L->!StringUtils.isBlank(L)).
						map(S->S.trim()).
						filter(S->!(S.equals("-") || S.equals(".") || S.startsWith("#"))).
						collect(Collectors.toSet());
					}
				}
			else
				{
				knownGeneIds = Collections.emptySet();
				}
			// open the sam file
			final VCFHeader header = iterin.getHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			final IntervalTreeMap<List<Intron>> intronMap = new IntervalTreeMap<>();
			
			try(CloseableIterator<UcscTranscript> iter = UcscTranscriptCodec.makeIterator(this.genPredURI)) {
					while(iter.hasNext()) {
						final UcscTranscript kg = iter.next();
						if(!kg.hasIntrons()) continue;
						if(dict.getSequence(kg.getContig())==null) continue;
						for(UcscTranscript.Intron kgi: kg.getIntrons()) {
							final Intron intron = new Intron(kgi);
							final Interval r = new Interval(intron);
							List<Intron> introns= intronMap.get(r);
							if(introns==null) {
								introns=new ArrayList<>();
								intronMap.put(r,introns);
								}
							introns.add(intron);
							}
						}
					}
			
			/** build vcf header */
			final VCFHeader header2=new VCFHeader(header); 
			header2.addMetaDataLine(new VCFFilterHeaderLine(ATT_FILTER_INTRON,"variant could be a deleted intron from a retrocopy"));
			header2.addMetaDataLine(new VCFFilterHeaderLine(ATT_FILTER_KNOWN,"variant could be a deleted intron from a known retrocopy"));
			header2.addMetaDataLine(new VCFInfoHeaderLine(ATT_RETROCOPY, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,"Identifiers for the retrocopies."));
			
			JVarkitVersion.getInstance().addMetaData(this, header2);
			
			out.writeHeader(header2);

			while(iterin.hasNext()){
				final VariantContext ctx= iterin.next();
				if(ctx.getStart()==ctx.getEnd()) {
					out.add(ctx);
					continue;
					}
				final String svType = ctx.getAttributeAsString(VCFConstants.SVTYPE, "");
				
				if(svType.equals("BND") || svType.equals("INS")) {
					out.add(ctx);
					continue;
					}
				
				boolean known_flag =false;
				final Set<String> retrocopy_identifiers= new TreeSet<>();
				for(final Intron intron: intronMap.getOverlapping(ctx).stream().
						flatMap(L->L.stream()).
						filter(I->isWithinDistance(I.getStart(),ctx.getStart())).
						filter(I->isWithinDistance(I.getEnd(),ctx.getEnd())).
						collect(Collectors.toList())
						) 
					{
					if(knownGeneIds.contains(intron.geneId)) {
						known_flag=true;
						}
					retrocopy_identifiers.add(VCFUtils.escapeInfoField(intron.geneId));
					retrocopy_identifiers.add(VCFUtils.escapeInfoField(intron.transcriptId));
					}					
				
				if(retrocopy_identifiers.isEmpty()) {
					out.add(ctx);
					continue;
					}
				
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(ATT_FILTER_INTRON);
				if(known_flag) vcb.filter(ATT_FILTER_KNOWN);
				vcb.attribute(ATT_RETROCOPY, new ArrayList<>(retrocopy_identifiers));
				
				out.add(vcb.make());
				}
						
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	public static void main(final String[] args) {
		new KnownRetroCopy().instanceMainWithExit(args);
	}
}
