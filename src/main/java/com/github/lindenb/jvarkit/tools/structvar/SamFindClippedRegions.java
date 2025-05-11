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
package com.github.lindenb.jvarkit.tools.structvar;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.function.BiFunction;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.structure.GtfReader;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
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
@Program(
	name="samfindclippedregions",
	description="Fins clipped position in one or more bam. ",
	keywords= {"sam","bam","clip","vcf"},
	creationDate="20140228",
	modificationDate="20220329",
	jvarkit_amalgamion = true,
	menu="CNV/SV"
	)
public class SamFindClippedRegions extends Launcher {
	private static final Logger LOG=Logger.of(SamFindClippedRegions.class);
	
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFai = null;
	@Parameter(names={"--bed","--regions-file"},description="restrict to this bed file. "+BedLineReader.OPT_DESC)
	private Path bedPath = null;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection =new WritingSortingCollection();
	


	@Parameter(names="-c",description="consider only reads with clip having length >= 'x'")
	private int min_clip_operator_length = 1;
	@Parameter(names={"--min-depth"},description="Ignore if Depth lower than 'x'")
	private int min_depth=10;
	@Parameter(names={"--min-clip-depth"},description="Ignore if number of clipped bases overlaping one POS is lower than 'x'")
	private int min_clip_depth =10;
	@Parameter(names={"--min-ratio"},description="Ignore genotypes where count(clip)/(count(clip)+DP) < x",splitter=NoSplitter.class,converter=FractionConverter.class)
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
	

	private static class Base
		{
		final int tid;
		final int pos;
		final int sample_idx;
		Base(final int sample_idx,final int tid,final int pos) {
			this.tid = tid;
			this.pos = pos;
			this.sample_idx = sample_idx;
			}
		
		
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
		int compare2(final Base o)  {
			int i = Integer.compare(this.tid, o.tid);
			if(i!=0) return i;
			return Integer.compare(this.pos, o.pos);
			}
		
		int compare1(final Base o)  {
			int i= compare2(o);
			if(i!=0) return i;
			return Integer.compare(sample_idx,o.sample_idx);
			}
		
		}
	
	private static class BaseCodec extends AbstractDataCodec<Base> {
		@Override
		public AbstractDataCodec<Base> clone() {
			return new BaseCodec();
			}
		@Override
		public Base decode(DataInputStream dis) throws IOException {
			int tid;
			try {
				tid = dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			int pos = dis.readInt();
			int sample_idx = dis.readInt();
			final Base b= new Base(sample_idx, tid, pos);
			b.noClip = dis.readInt();
			b.leftClip = dis.readInt();
			b.rightClip = dis.readInt();
			b.del = dis.readInt();
			b.noClip_sum_mapq = dis.readDouble();
			return b;
			}
		@Override
		public void encode(DataOutputStream dos, Base o) throws IOException {
			dos.writeInt(o.tid);
			dos.writeInt(o.pos);
			dos.writeInt(o.sample_idx);
			
			dos.writeInt(o.noClip);
			dos.writeInt(o.leftClip);
			dos.writeInt(o.rightClip);
			dos.writeInt(o.del);
			dos.writeDouble(o.noClip_sum_mapq);
			}
	}
	
	
	@Override
	public int doWork(final List<String> args) {
		final int bad_mapq = 30;
		try {
			if(this.min_clip_depth >this.min_depth) {
				LOG.error("this.min_clip_depth>this.min_depth");
				return -1;
			}
			if(this.fraction<0 || this.fraction>1.0) {
				LOG.error("bad ratio: "+fraction);
				return -1;
			}
			
			final SAMSequenceDictionary dict =  SequenceDictionaryUtils.extractRequired(this.referenceFai);
			final List<Path> inputs = IOUtils.unrollPaths(args);
			if(inputs.isEmpty()) {
				LOG.error("input is missing");
				return -1;
			}
			
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
			
			
			final List<String> sample_list = new ArrayList<>();
			
			
			final QueryInterval[] queryIntervals;
			if(this.bedPath==null) {
				queryIntervals = null;
				} else {
				try(BedLineReader br = new BedLineReader(this.bedPath)){
					br.setValidationStringency(ValidationStringency.LENIENT).
						setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
					queryIntervals = br.optimizeIntervals(dict);
					}
				}
			
			SortingCollection<Base> sortingCollection =
					SortingCollection.newInstance(
							Base.class,
							new BaseCodec(),
							(A,B)->A.compare1(B),
							writingSortingCollection.getMaxRecordsInRam(),
							writingSortingCollection.getTmpPath()
							);
			sortingCollection.setDestructiveIteration(true);
			
			final Predicate<Base> acceptBase = B->{
				return B.clip()>0;
				};
			
			final SamReaderFactory srf = super.createSamReaderFactory().referenceSequence(this.referenceFai);
			for(final Path path: inputs) {
				try(SamReader sr = srf.open(path)) {
					final SAMFileHeader header = sr.getFileHeader();
					SequenceUtil.assertSequenceDictionariesEqual(dict, SequenceDictionaryUtils.extractRequired(header));

					final String sample_name = header.getReadGroups().
						stream().
						map(RG->RG.getSample()).
						filter(S->!StringUtils.isBlank(S)).
						findFirst().
						orElse(IOUtils.getFilenameWithoutCommonSuffixes(path))
						;
					if(sample_list.contains(sample_name)) {
						LOG.error("duplicate sample "+sample_name+" in "+path);
						return -1;
						}
					final int sample_idx = sample_list.size();
					sample_list.add(sample_name);
					
					int prev_tid = -1;
					final SortedMap<Integer,Base> pos2base = new TreeMap<>();
					/* get base in pos2base, create it if needed */
					final BiFunction<Integer,Integer, Base> baseAt = (TID,POS)->{
						Base b = pos2base.get(POS);
						if(b==null) {
							b = new Base(sample_idx,TID,POS);
							pos2base.put(POS, b);
							}
						return b;
						};
					
					
					try(CloseableIterator<SAMRecord> iter= queryIntervals==null?sr.iterator():sr.query(queryIntervals, false)) {
						for(;;) {
							final SAMRecord rec=(iter.hasNext()?iter.next():null);
							if(rec!=null && !SAMRecordDefaultFilter.accept(rec,this.min_mapq)) continue;
							
							if(rec==null || prev_tid!=rec.getReferenceIndex().intValue())
								{
								for(final Integer pos: pos2base.keySet()) {
									final Base b = pos2base.get(pos);
									if(acceptBase.test(b)) sortingCollection.add(b);
									}
								if(rec==null) break;
								pos2base.clear();
								prev_tid = rec.getReferenceIndex().intValue();
								}
							/* pop back previous bases */
							for(Iterator<Integer> rpos = pos2base.keySet().iterator();
									rpos.hasNext();
								) {
						    	final Integer pos = rpos.next();
						    	if(pos.intValue() + this.max_clip_length >= rec.getUnclippedStart()) break;
						    	final Base b = pos2base.get(pos);
						    	if(acceptBase.test(b)) sortingCollection.add(b);
						    	rpos.remove();
						    	}
							
							final Cigar cigar = rec.getCigar();
							int refPos= rec.getAlignmentStart();
							for(final CigarElement ce: cigar.getCigarElements()) {
								final CigarOperator op = ce.getOperator();
								if(op.consumesReferenceBases()) {
									if(op.consumesReadBases()) {
										for(int x=0;x< ce.getLength();++x) {
											final Base gt=baseAt.apply(prev_tid,refPos+x);
											gt.noClip++;
											gt.noClip_sum_mapq += rec.getMappingQuality();
											}
										}
									else if(op.equals(CigarOperator.D) || op.equals(CigarOperator.N))
										{
										baseAt.apply(prev_tid,refPos).del++;
										baseAt.apply(prev_tid,refPos + ce.getLength() - 1).del++;
										}
									refPos += ce.getLength();
									}
								}
							
							CigarElement ce = cigar.getFirstCigarElement();
							if(ce!=null && ce.getOperator().isClipping() && ce.getLength()>=this.min_clip_operator_length) {
								baseAt.apply(prev_tid,rec.getStart()-1).leftClip++;
								}
							ce = cigar.getLastCigarElement();
							if(ce!=null && ce.getOperator().isClipping() && ce.getLength()>=this.min_clip_operator_length) {
								baseAt.apply(prev_tid,rec.getEnd()+1).rightClip++;
								}
							
							}
						}
					/* write last bases */
					for(final Integer pos: pos2base.keySet()) {
						final Base b = pos2base.get(pos);
						if(acceptBase.test(b)) sortingCollection.add(b);
						}
					} // end open reader
				}//end loop sam files
		
			sortingCollection.doneAdding();
			

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

			
			final VCFHeader vcfHeader=new VCFHeader(vcfHeaderLines,sample_list);
			vcfHeader.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
			
			this.writingVcfConfig.dictionary(dict);
			try(VariantContextWriter w = this.writingVcfConfig.open(this.outputFile)) {
				w.writeHeader(vcfHeader);
				try(CloseableIterator<Base> r1=sortingCollection.iterator()) {
					try(EqualRangeIterator<Base> r2=new EqualRangeIterator<>(r1,(A,B)->A.compare2(B))) {
						while(r2.hasNext()) {
							final List<Base> array = r2.next();
							final Base first = array.get(0);
							if(first.pos<1) continue;
							
							//no clip
							if(array.stream().mapToInt(G->G.clip()).sum()==0) continue;
							
							if(array.stream().allMatch(G->G.clip() < min_clip_depth)) continue;
							if(array.stream().allMatch(G->G.dp() < min_depth)) continue;
							
							
							if(array.stream().allMatch(G->G.ratio() < fraction)) continue;
							final VariantContextBuilder vcb=new VariantContextBuilder();
							final String chrom = dict.getSequence(first.tid).getSequenceName();
							vcb.chr(chrom);
							vcb.start(first.pos);
							vcb.stop(first.pos);
							vcb.alleles(Arrays.asList(reference_allele,alt_allele));
							vcb.attribute(VCFConstants.DEPTH_KEY,array.stream().mapToInt(G->G.dp()).sum());
							
							

							
							final List<Genotype> genotypes = new ArrayList<>(array.size());
							int AC=0;
							int AN=0;
							int max_clip=1;
							double sum_mapq=0.0;
							int count_mapq = 0;
							
							for(final Base gt : array) {
								final GenotypeBuilder gb = new GenotypeBuilder(sample_list.get(gt.sample_idx));
								
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
							
							if(gtfPath!=null) {
								final Locatable bounds1 = new SimpleInterval(chrom, Math.max(1, first.pos-max_intron_distance), first.pos+max_intron_distance);
								intronMap.getOverlapping(bounds1).stream().
									filter(I->
											Math.abs(I.getStart()-first.pos)<= this.max_intron_distance ||
											Math.abs(I.getEnd()-first.pos)<= this.max_intron_distance
											).
									map(I->I.getName()).
									findFirst().
									ifPresent(transcript_id->{
										vcb.attribute(infoRetrogene.getID(),transcript_id);
										vcb.filter(filterRetrogene.getID());
										});
										;						
								}
							
							vcb.log10PError(max_clip/-10.0);
							vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, AC);
							vcb.attribute(VCFConstants.ALLELE_NUMBER_KEY, AN);
							if(AN>0) vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, AC/(float)AN);
							vcb.genotypes(genotypes);
							w.add(vcb.make());
							}//end while r2
						}//end r2
					}//end r1
				}//end writer
			sortingCollection.cleanup();
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	public static void main(final String[] args)
		{
		new SamFindClippedRegions().instanceMainWithExit(args);
		}
	}
