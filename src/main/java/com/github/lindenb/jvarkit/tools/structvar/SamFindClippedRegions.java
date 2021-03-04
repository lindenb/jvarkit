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
package com.github.lindenb.jvarkit.tools.structvar;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.jcommander.MultiBamLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

samfindclippedregions find 'blunt' regions where the reads are clipped. It can be used to find structural variations in exomes/genomes data.

input is a set of indexed BAM/CRAM files or a list with the '.list' suffix containing the path to the bam.

output is a VCF file


### Example

```
$ java -jar dist/samfindclippedregions.jar --min-depth 10 --min-ratio 0.2 src/test/resources/S*.bam


##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=CL,Number=1,Type=Integer,Description="Left Clip">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RL,Number=1,Type=Integer,Description="Right Clip">
##FORMAT=<ID=TL,Number=1,Type=Integer,Description="Total Clip">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
(...)
##samfindclippedregions.meta=compilation:20191009163450 githash:b8d60cab htsjdk:2.20.1 date:20191009163608 cmd:--min-depth 10 --min-ratio 0.2 src/test/resources/S1.bam src/test/resources/S2.bam src/test/resources/S3.bam src/test/resources/S4.bam src/test/resources/S5.bam
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5
RF01	996	.	N	<CLIP>	.	.	AC=1;AF=0.1;AN=10;DP=30	GT:AD:CL:DP:RL:TL	0/0:2,0:0:2:0:0	0/0:4,0:0:4:0:00/0:4,0:0:4:0:0	0/0:15,0:0:15:0:0	0/1:4,1:0:5:1:1
(...)
```


### Screenshot

https://twitter.com/yokofakun/status/1194921855875977216

![https://twitter.com/yokofakun/status/1194921855875977216](https://pbs.twimg.com/media/EJU3F9hWoAACgsd?format=png&name=large)

END_DOC

 */
@Program(name="samfindclippedregions",
	description="Fins clipped position in one or more bam. ",
	keywords= {"sam","bam","clip","vcf"},
	creationDate="20140228",
	modificationDate="20210304"
	)
public class SamFindClippedRegions extends MultiBamLauncher
	{
	private static final Logger LOG=Logger.build(SamFindClippedRegions.class).make();
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names="-c",description="consider only clip having length >= 'x'")
	private int min_clip_operator_length = 1;
	@Parameter(names={"--groupby","--partition"},description="Group Reads by. "+SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition=SAMRecordPartition.sample;
	@Parameter(names={"--min-depth"},description="Ignore if Depth lower than 'x'")
	private int min_depth=10;
	@Parameter(names={"--min-clip-depth"},description="Ignore if number of clipped bases lower than 'x'")
	private int min_clip_depth =10;
	@Parameter(names={"--min-ratio"},description="Ignore genotypes where count(clip)/(count(clip)+DP) < x")
	private double fraction = 0.1;
	@Parameter(names={"--gtf"},description="Optional gtf file. Will be used to set a warning if the junction could be a junction exon-exon of a retrogene. "+GtfReader.OPT_DESC)
	private Path gtfPath = null;
	@Parameter(names={"--intron-distance"},description="when gtf is specified: max distance between breakend and the intron bound")
	private int max_intron_distance=3;
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int min_mapq= 1;

	@ParametersDelegate
	private WritingVariantsDelegate writingVcfConfig = new WritingVariantsDelegate();
	
	private final int max_clip_length=1000;
	
	private static class Gt
		{
		int noClip=0;
		int leftClip = 0;
		int rightClip = 0;
		int del = 0;
		double noClip_sum_mapq = 0;
		int clip() { return leftClip+rightClip;}
		int dp() { return noClip+clip();}
		double ratio() {return clip()/(double)dp();}
		int noClipMapq() {
			return noClip==0 ? 0: (int) (this.noClip_sum_mapq /this.noClip);
			}
		}

	private static class Base
		{
		final int pos;
		final Map<String,Gt> sample2gt;
		Base(final int pos,final Set<String> samples) {
			this.pos = pos;
			sample2gt = new HashMap<>(samples.size());
			for(final String sn:samples) sample2gt.put(sn, new Gt());
			}
		Gt getGt(final String sn) {
			return Objects.requireNonNull(this.sample2gt.get(sn),sn);
			}
		}
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int beforeSam() {
		if(this.min_clip_depth >this.min_depth) {
			LOG.error("this.min_clip_depth>this.min_depth");
			return -1;
		}
		if(this.fraction<0 || this.fraction>1.0) {
			LOG.error("bad ratio: "+fraction);
			return -1;
		}
		return super.beforeSam();
		}
	
	
	@Override
	protected int processInput(final SAMFileHeader header, final CloseableIterator<SAMRecord> iter) {
		final Set<String> samples = header.getReadGroups().
			stream().
			map(this.partition).
			filter(S->!StringUtils.isBlank(S)).
			collect(Collectors.toCollection(TreeSet::new));
	
		if(samples.isEmpty()) {
			LOG.error("no sample in read group was defined.");
			return -1;
		}
		
		final SAMSequenceDictionary dict =  SequenceDictionaryUtils.extractRequired(header);
		
		final int bad_mapq = 30;
		
		//SAMFileWriter w=null;
		try
			{			
			
			final IntervalTreeMap<Interval> intronMap = new IntervalTreeMap<>();
			if(this.gtfPath!=null) {
				try(GtfReader gtfReader= new GtfReader(this.gtfPath)) {
					gtfReader.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					gtfReader.getAllGenes().
						stream().
						flatMap(G->G.getTranscripts().stream()).
						filter(T->T.hasIntron()).
						flatMap(T->T.getIntrons().stream()).
						map(I->new Interval(I.getContig(),I.getStart(),I.getEnd(),I.isNegativeStrand(),I.getTranscript().getId())).
						forEach(I->intronMap.put(I,I));
				}
			}
			
			/* build VCF header */
			final Allele reference_allele= Allele.create("N",true);			
			final Allele alt_allele = Allele.create("<CLIP>",false);			
		
			
			final Set<VCFHeaderLine> vcfHeaderLines=new HashSet<>();
			
			VCFStandardHeaderLines.addStandardFormatLines(vcfHeaderLines, true, 
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.DEPTH_KEY,
					VCFConstants.GENOTYPE_ALLELE_DEPTHS
					);
			
			VCFStandardHeaderLines.addStandardInfoLines(vcfHeaderLines, true, 
					VCFConstants.DEPTH_KEY,
					VCFConstants.ALLELE_COUNT_KEY,
					VCFConstants.ALLELE_NUMBER_KEY,
					VCFConstants.ALLELE_FREQUENCY_KEY
					)
					;

			final VCFFormatHeaderLine leftClip = new VCFFormatHeaderLine("CL", 1,VCFHeaderLineType.Integer,"Left Clip");
			vcfHeaderLines.add(leftClip);
			final VCFFormatHeaderLine rightClip = new VCFFormatHeaderLine("RL", 1,VCFHeaderLineType.Integer,"Right Clip");
			vcfHeaderLines.add(rightClip);
			final VCFFormatHeaderLine totalCip = new VCFFormatHeaderLine("TL", 1,VCFHeaderLineType.Integer,"Total Clip");
			vcfHeaderLines.add(totalCip);
			final VCFFormatHeaderLine totalDel = new VCFFormatHeaderLine("DL", 1,VCFHeaderLineType.Integer,"Total Deletions");
			vcfHeaderLines.add(totalDel);
			final VCFFormatHeaderLine noClipMAPQ = new VCFFormatHeaderLine("MQ", 1,VCFHeaderLineType.Integer,"Average MAPQ for reads without clip at this position.");
			vcfHeaderLines.add(noClipMAPQ);

			final VCFInfoHeaderLine averageMAPQ = new VCFInfoHeaderLine("AVG_MAPQ", 1,VCFHeaderLineType.Integer,"Average MAPQ for called genotypes");
			vcfHeaderLines.add(averageMAPQ);

			final VCFInfoHeaderLine infoRetrogene = new VCFInfoHeaderLine("RETROGENE", 1,VCFHeaderLineType.String,"transcript name for Possible retrogene.");
			vcfHeaderLines.add(infoRetrogene);
			final VCFFilterHeaderLine filterRetrogene = new VCFFilterHeaderLine("POSSIBLE_RETROGENE","Junction is a possible Retrogene.");
			vcfHeaderLines.add(filterRetrogene);
			final VCFFilterHeaderLine filterlowMapq = new VCFFilterHeaderLine("LOW_MAPQ","Low average mapq (< "+bad_mapq+")");
			vcfHeaderLines.add(filterlowMapq);

			
			final VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,samples);
			vcfHeader.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
			
			this.writingVcfConfig.dictionary(dict);
				try(VariantContextWriter w = this.writingVcfConfig.open(this.outputFile)) {
				
			
				w.writeHeader(vcfHeader);
				
				@SuppressWarnings("resource")
				final VariantContextWriter finalVariantContextWriter = w;
				
				/** dump a BASe into the VCF */
				final BiConsumer<String,Base> baseConsumer = (CTG,B)->{
					if(B.pos<1) return;
				
					//no clip
					if(B.sample2gt.values().stream().mapToInt(G->G.clip()).sum()==0) return;
					
					if(B.sample2gt.values().stream().allMatch(G->G.clip() < min_clip_depth)) return;
					if(B.sample2gt.values().stream().allMatch(G->G.dp() < min_depth)) return;
					
					
					if(B.sample2gt.values().stream().allMatch(G->G.ratio() < fraction)) return;
					final VariantContextBuilder vcb=new VariantContextBuilder();
					vcb.chr(CTG);
					vcb.start(B.pos);
					vcb.stop(B.pos);
					vcb.alleles(Arrays.asList(reference_allele,alt_allele));
					vcb.attribute(VCFConstants.DEPTH_KEY,B.sample2gt.values().stream().mapToInt(G->G.dp()).sum());
					
					/* if gtf was specified, find intron which ends are near this pos */
					if(gtfPath!=null) {
						final Locatable bounds1 = new SimpleInterval(CTG, Math.max(1, B.pos-max_intron_distance), B.pos+max_intron_distance);
						intronMap.getOverlapping(bounds1).stream().
							filter(I->
									Math.abs(I.getStart()-B.pos)<= this.max_intron_distance ||
									Math.abs(I.getEnd()-B.pos)<= this.max_intron_distance
									).
							map(I->I.getName()).
							findFirst().
							ifPresent(transcript_id->{
								vcb.attribute(infoRetrogene.getID(),transcript_id);
								vcb.filter(filterRetrogene.getID());
								});
								;						
						}
					
					
					final List<Genotype> genotypes = new ArrayList<>(B.sample2gt.size());
					int AC=0;
					int AN=0;
					int max_clip=1;
					double sum_mapq=0.0;
					int count_mapq = 0;
					
					for(final String sn:B.sample2gt.keySet()) {
						final Gt gt = B.sample2gt.get(sn);
						final GenotypeBuilder gb = new GenotypeBuilder(sn);
						
						if(gt.clip()==0 && gt.noClip==0)
							{
							gb.alleles(Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));
							}
						else if(gt.noClip==0) {
							gb.alleles(Arrays.asList(alt_allele,alt_allele));
							AC+=2;
							sum_mapq += gt.noClipMapq();
							count_mapq++;
							AN+=2;
							}
						else if(gt.clip()==0) {
							gb.alleles(Arrays.asList(reference_allele,reference_allele));
							AN+=2;
							}
						else{
							gb.alleles(Arrays.asList(reference_allele,alt_allele));
							AC++;
							sum_mapq += gt.noClipMapq();
							count_mapq++;
							AN+=2;
							}
						
						gb.DP(gt.dp());
						gb.attribute(leftClip.getID(), gt.leftClip);
						gb.attribute(rightClip.getID(), gt.rightClip);
						gb.attribute(totalCip.getID(), gt.clip());
						gb.attribute(totalDel.getID(), gt.del);
						gb.attribute(noClipMAPQ.getID(), gt.noClipMapq());
						gb.AD(new int[] {gt.noClip, gt.clip()});
						
						genotypes.add(gb.make());
						
						max_clip = Math.max(max_clip, gt.clip());
						}
					if(count_mapq>0) {
						final int avg_mapq = (int)(sum_mapq/count_mapq);
						vcb.attribute(averageMAPQ.getID(),avg_mapq);
						if(avg_mapq < bad_mapq ) vcb.filter(filterlowMapq.getID());
						}
					vcb.log10PError(max_clip/-10.0);
					vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, AC);
					vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY, AN);
					if(AN>0) vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, AC/(float)AN);
					vcb.genotypes(genotypes);
					finalVariantContextWriter.add(vcb.make());
					};
					
					
				final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
						newInstance().
						dictionary(dict).
						logger(LOG).
						build();
				
				String prevContig = null;
				final SortedMap<Integer,Base> pos2base = new TreeMap<>();
				
				/* get base in pos2base, create it if needed */
				final Function<Integer, Base> baseAt = POS->{
					Base b = pos2base.get(POS);
					if(b==null) {
						b = new Base(POS,samples);
						pos2base.put(POS, b);
						}
					return b;
					};
				
				for(;;) {
					final SAMRecord rec=(iter.hasNext()?progress.apply(iter.next()):null);
					if(rec!=null && !SAMRecordDefaultFilter.accept(rec,this.min_mapq)) continue;
					
					if(rec==null || !rec.getContig().equals(prevContig))
						{
						for(final Integer pos: pos2base.keySet()) {
							baseConsumer.accept(prevContig,pos2base.get(pos));
							}
						if(rec==null) break;
						pos2base.clear();
						prevContig = rec.getContig();
						}
					
					for(Iterator<Integer> rpos = pos2base.keySet().iterator();
							rpos.hasNext();
						) {
				    	final Integer pos = rpos.next();
				    	if(pos.intValue() + this.max_clip_length >= rec.getUnclippedStart()) break;
				    	baseConsumer.accept(prevContig,pos2base.get(pos));
				    	rpos.remove();
				    	}
	
					
					
					final String rg = this.partition.getPartion(rec);
					if(StringUtils.isBlank(rg)) continue;
					
					for(final AlignmentBlock ab:rec.getAlignmentBlocks()) {
						for(int n=0;n< ab.getLength();++n) {
							
							}
						}
					
					final Cigar cigar = rec.getCigar();
					int refPos= rec.getAlignmentStart();
					for(final CigarElement ce: cigar.getCigarElements()) {
						final CigarOperator op = ce.getOperator();
						if(op.consumesReferenceBases()) {
							if(op.consumesReadBases()) {
								for(int x=0;x< ce.getLength();++x) {
									final Gt gt=baseAt.apply(refPos+x).getGt(rg);
									gt.noClip++;
									gt.noClip_sum_mapq += rec.getMappingQuality();
									}
								}
							else if(op.equals(CigarOperator.D) || op.equals(CigarOperator.N))
								{
								baseAt.apply(refPos).getGt(rg).del++;
								baseAt.apply(refPos + ce.getLength() - 1).getGt(rg).del++;
								}
							refPos += ce.getLength();
							}
						}
					
					CigarElement ce = cigar.getFirstCigarElement();
					if(ce!=null && ce.getOperator().isClipping() && ce.getLength()>=this.min_clip_operator_length) {
						baseAt.apply(rec.getStart()-1).getGt(rg).leftClip++;
						}
					ce = cigar.getLastCigarElement();
					if(ce!=null && ce.getOperator().isClipping() && ce.getLength()>=this.min_clip_operator_length) {
						baseAt.apply(rec.getEnd()+1).getGt(rg).rightClip++;
						}
					
					}
				
				}// end of vcf writer
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args)
		{
		new SamFindClippedRegions().instanceMainWithExit(args);
		}
	}
