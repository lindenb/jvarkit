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
package com.github.lindenb.jvarkit.tools.retrocopy;

import java.io.BufferedReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.Gene;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.bio.structure.Intron;
import com.github.lindenb.jvarkit.util.bio.structure.Transcript;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Example

```
java -jar dist/gtfretrocopy.jar --gtf transcript.gtf.gz input.vcf.gz > retrocopies.vcf
```

END_DOC

*/
@Program(name="gtfretrocopy",
description="Scan retrocopies by comparing the gtf/intron and the deletions in a VCF",
keywords={"gtf","retrocopy","deletion"},
creationDate="20190813",
modificationDate="20191104"
)
public class GtfRetroCopy extends Launcher
	{
	private static final Logger LOG = Logger.build(GtfRetroCopy.class).make();
	private enum IdKey {transcript_id,gene_id,gene_name};
	
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-gtf","--gtf"},description=GtfReader.OPT_DESC,required=true)
	private Path gtfPath = null;
	@Parameter(names={"-d","--distance"},description="max distance between an intron and the deletion found in the VCF")
	private int distance = 10;
	@Parameter(names={"--mic","--min-intron-count"},description="Min intron count.",hidden=true)
	private int min_intron_count = 0;
	@Parameter(names={"--all"},description="all introns must be found")
	private boolean only_all_introns = false;
	@Parameter(names={"-k","--known"},description="Gene-ID of known retrogenes. One per line. A source could be : http://retrogenedb.amu.edu.pl/static/download/")
	private Path knownPath = null;
	@Parameter(names={"--id","-id"},description="Which key should I use for the column ID. The idea is to use the gene name to get the uniq entities per vcf.")
	private IdKey idKey = IdKey.transcript_id;

	@ParametersDelegate
	private WritingVariantsDelegate writingVcf = new WritingVariantsDelegate();

	
	
	
	private static final String ENSEMBL_TRANSCRIPT_ATTS[]=new String[] {"gene_id","gene_version","transcript_id","transcript_version","gene_name","gene_source","gene_biotype","transcript_name","transcript_source","transcript_biotype","tag","ccds_id","havana_transcript","havana_transcript_version","tag"};
	private final static String ATT_INTRONS_COUNT="COUNT_INTRONS";
	private final static String ATT_INTRONS_CANDIDATE_COUNT="MATCHING_INTRONS_COUNT";
	private final static String ATT_INTRONS_CANDIDATE_FRACTION="MATCHING_INTRONS_FRACTION";
	private final static String ATT_NOT_ALL_INTRONS="NOT_ALL_INTRONS";
	private final static String ATT_INTRONS_BOUNDS="INTRONS_BOUNDS";
	private final static String ATT_INTRONS_SIZES="INTRONS_SIZES";
	private final static String ATT_KNOWN="KNOWN_RETROGENE";

	
	private boolean isWithinDistance(int a, int b) {
		return Math.abs(a-b) < this.distance;
	}
	
	@Override
	public int doWork(final List<String> args) {
		VCFReader vcfFileReader = null;
		VariantContextWriter vcw0=null;
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
			final String input = oneAndOnlyOneFile(args);
			vcfFileReader = VCFReaderFactory.makeDefault().open(Paths.get(input), true);
			final VCFHeader header = vcfFileReader.getHeader();
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			
			final Comparator<String> contigCmp =
					dict == null?
					(A,B)->A.compareTo(B):
					new ContigDictComparator(dict)
					;
			final Comparator<Gene> geneCmp = (A,B)->{
				int i= contigCmp.compare(A.getContig(),B.getContig());
				if(i!=0) return i;
				i =  Integer.compare(A.getStart(),B.getStart());
				if(i!=0) return i;
				return Integer.compare(A.getEnd(),B.getEnd());
			};
			
			final GtfReader gtfReader = new GtfReader(this.gtfPath);
			if(dict!=null && !dict.isEmpty()) {
				this.writingVcf.dictionary(dict);
				gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
			}
			final List<Gene> genes = gtfReader.getAllGenes().
					stream().
					filter(G->G.getTranscripts().stream().count()>0L).
					filter(G->G.getTranscripts().stream().anyMatch(T->T.getIntronCount()>=this.min_intron_count)).
					sorted(geneCmp).
					collect(Collectors.toList())
					;
			gtfReader.close();
			
			/** build vcf header */
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_QUALITY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_ALLELE_DEPTHS,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_NUMBER_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_COUNT_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.ALLELE_FREQUENCY_KEY,true));
			metaData.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.END_KEY,true));
			metaData.add(new VCFInfoHeaderLine(VCFConstants.SVTYPE, 1, VCFHeaderLineType.String,"Variation type"));
			metaData.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"Variation Length"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_BOUNDS, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,"Introns boundaries"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_SIZES, VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.Integer,"Introns sizes"));
			metaData.add(new VCFFilterHeaderLine(ATT_NOT_ALL_INTRONS,"Not all introns were found retrocopied"));

			for(final String att:ENSEMBL_TRANSCRIPT_ATTS)
				{
				metaData.add(new VCFInfoHeaderLine(att, 1, VCFHeaderLineType.String,"Value for the attribute '"+att+"' in the gtf"));
				}
			
			
			//metaData.add(new VCFFormatHeaderLine(ATT_COUNT_SUPPORTING_READS, 2,VCFHeaderLineType.Integer,"Count supporting reads [intron-left/intron-right]"));
			//metaData.add(new VCFInfoHeaderLine(ATT_RETRO_DESC, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,
			//		"Retrocopy attributes: transcript-id|strand|exon-left|exon-left-bases|exon-right-bases|exon-right"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_COUNT, 1, VCFHeaderLineType.Integer,"Number of introns for the Transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_CANDIDATE_COUNT, 1, VCFHeaderLineType.Integer,"Number of introns found retrocopied for the transcript"));
			metaData.add(new VCFInfoHeaderLine(ATT_INTRONS_CANDIDATE_FRACTION, 1, VCFHeaderLineType.Float,"Fraction of introns found retrocopied for the transcript"));
			metaData.add(new VCFFilterHeaderLine(ATT_NOT_ALL_INTRONS,"Not all introns were found retrocopied"));
			metaData.add(new VCFFilterHeaderLine(ATT_KNOWN,"Known RetroGenes. "+(this.knownPath==null?"":" Source: "+this.knownPath)));
			
			
			
			
			
			final VCFHeader header2=new VCFHeader(header); 
			metaData.stream().forEach(H->header2.addMetaDataLine(H));
			
			JVarkitVersion.getInstance().addMetaData(this, header2);
			
			final Allele ref= Allele.create((byte)'N', true);
			final Allele alt= Allele.create("<RETROCOPY>", false);
			
			
			/* open vcf for writing*/
			vcw0= this.writingVcf.open(this.outputFile);
			vcw0.writeHeader(header2);

			
			final ProgressFactory.Watcher<Gene> progress = ProgressFactory.newInstance().logger(LOG).dictionary(dict).build();
			for(final Gene gene : genes) {
				progress.apply(gene);
				final List<VariantContext> variants = new ArrayList<>();
				final CloseableIterator<VariantContext> iter2 = vcfFileReader.query(gene.getContig(),gene.getStart(),gene.getEnd());
				while(iter2.hasNext()) {
					final VariantContext ctx = iter2.next();
					if(ctx.getStart()==ctx.getEnd()) continue;//SNV
					StructuralVariantType svType = ctx.getStructuralVariantType();
					if(StructuralVariantType.BND.equals(svType)) continue;
					if(StructuralVariantType.INS.equals(svType)) continue;
					variants.add(ctx);
					}
				iter2.close();
				if(variants.isEmpty()) continue;
				
				for(final Transcript transcript:gene.getTranscripts()) {
					if(!transcript.hasIntron()) continue;
					if(transcript.getIntronCount() < this.min_intron_count) continue;
					final Counter<String> samples = new Counter<>();
					for(final Intron intron: transcript.getIntrons()) {
						for(final VariantContext ctx:variants) {
							if(!isWithinDistance(intron.getStart(),ctx.getStart())) continue;
							if(!isWithinDistance(intron.getEnd(),ctx.getEnd())) continue;
							if(ctx.hasGenotypes()) {
								for(final Genotype g:ctx.getGenotypes()) {
									if(g.isNoCall() || g.isHomRef()) continue;
									samples.incr(g.getSampleName());
									}
								}
							else
								{
								samples.incr("*");
								}
							} // end iter2
						}// end intron
					final long max_count = samples.stream().mapToLong(E->E.getValue()).max().orElse(0L);
					if(max_count==0) continue;
					if(this.only_all_introns && max_count!=transcript.getIntronCount()) continue;
					
					// ok good candidate
					final VariantContextBuilder vcb = new VariantContextBuilder();
					vcb.chr(transcript.getContig());
					vcb.start(transcript.getStart());
					vcb.stop(transcript.getEnd());
					switch(this.idKey) {
						case gene_name:
							final String gn= transcript.getGene().getGeneName();
							vcb.id( StringUtils.isBlank(gn)?transcript.getId():gn);
							break;
						case gene_id:vcb.id(transcript.getGene().getId());break;
						case transcript_id : vcb.id(transcript.getId());break;
						default: throw new IllegalStateException();
						}
					
					final List<Allele> alleles = Arrays.asList(ref,alt);

					//vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY,2);
					//vcb.attribute(VCFConstants.ALLELE_COUNT_KEY,1);
					//vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY,0.5);
					vcb.attribute(VCFConstants.SVTYPE,"DEL");
					vcb.attribute(VCFConstants.END_KEY,transcript.getEnd());
					vcb.attribute("SVLEN",transcript.getLengthOnReference());
					vcb.attribute(ATT_INTRONS_BOUNDS, transcript.getIntrons().stream().map(S->""+S.getStart()+"-"+S.getEnd()).collect(Collectors.toList()));
					vcb.attribute(ATT_INTRONS_SIZES, transcript.getIntrons().stream().mapToInt(S->S.getLengthOnReference()).toArray());
					
					
					for(final String att:ENSEMBL_TRANSCRIPT_ATTS) {
						final String v= transcript.getProperties().get(att);
						if(StringUtils.isBlank(v)) continue;
						vcb.attribute(att, v);
						}
					
					vcb.alleles(alleles);
					boolean pass_filter=true;
					// introns sequences
					vcb.attribute(ATT_INTRONS_CANDIDATE_COUNT,max_count);
					vcb.attribute(ATT_INTRONS_COUNT,transcript.getIntronCount());
					vcb.attribute(ATT_INTRONS_CANDIDATE_FRACTION,max_count/(float)transcript.getIntronCount());
					if(transcript.getIntronCount()!=max_count) {
						vcb.filter(ATT_NOT_ALL_INTRONS);
						pass_filter=false;
						}
					if(knownGeneIds.contains(transcript.getGene().getId())) {
						vcb.filter(ATT_KNOWN);
						pass_filter=false;
						}
					
					if(header.hasGenotypingData()) {
						final List<Genotype> genotypes = new ArrayList<>();
						for(final String sn: header.getSampleNamesInOrder()) {
							final List<Allele> gtalleles;
							if(samples.count(sn)==0L)
								{
								gtalleles=Arrays.asList(ref,ref);
								}
							else
								{
								gtalleles=Arrays.asList(ref,alt);
								}	
							final GenotypeBuilder gb = new GenotypeBuilder(sn,gtalleles);
							genotypes.add(gb.make());
							}
						vcb.genotypes(genotypes);
						}					
					if(pass_filter) vcb.passFilters();
					vcw0.add(vcb.make());
				}
			}
						
			progress.close();
			vcw0.close();
			vcfFileReader.close();
			vcfFileReader=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfFileReader);
			CloserUtil.close(vcw0);
			}
		}
	
	public static void main(final String[] args) {
		new GtfRetroCopy().instanceMainWithExit(args);
	}
}
