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

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.TreeSet;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.ArchiveFactory;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.math.DiscreteMedian;
import com.github.lindenb.jvarkit.math.stats.Percentile;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.tools.structvar.indexcov.IndexCovUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
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
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import htsjdk.tribble.readers.TabixReader;
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
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**
BEGIN_DOC

## Example

```
 java -jar dist/coveragematrix.jar -R src/test/resources/rotavirus_rf.fa --exclude gaps.tsv.gz src/test/resources/S*.bam
```


END_DOC 
 */
@Program(
	name="coveragematrix",
	description="generate a VCF file from bam coverage",
	keywords={"cnv","bam","depth","coverage"},
	creationDate="20200618",
	modificationDate="20200618",
	generate_doc=false
	)
public class CoverageMatrix extends Launcher {
	private static final Logger LOG = Logger.build( CoverageMatrix.class).make();
	@Parameter(names={"-o","--output"},description=ArchiveFactory.OPT_DESC)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path refPath = null;
	@Parameter(names={"--mapq"},description = "min mapping quality")
	private int min_mapq=1;
	@Parameter(names={"--max-depth"},description = "ignore position if depth > 'x'")
	private int max_depth=500;
	@Parameter(names = {"--black","--exclude"}, description = "Optional. BED Tabix indexed black-listed region")
	private Path blackListedPath=null;
	@Parameter(names = {"--bin","--bin-size"}, description = "Bin size. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int bin_size = 1_000;
	@Parameter(names = {"--chrom","--contig"}, description = "Restrict to that contig.")
	private String restrictContig = null;
	@Parameter(names = {"--treshold"}, description = IndexCovUtils.TRESHOLD_OPT_DESC)
	private double indexCovTreshold = IndexCovUtils.DEFAULT_TRESHOLD;
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection= new WritingSortingCollection();
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();


	

	private static class CovItem {
		int sample_idx;
		int pos;
		float depth;
		float stddev;
		int compare0(CovItem item) {
			int i = Integer.compare(this.pos, item.pos);
			if(i!=0) return i;
			return  Integer.compare(this.sample_idx, item.sample_idx);
			}
		}

	private static class CovItemCodec extends AbstractDataCodec<CovItem> {
		@Override
		public CovItem decode(final DataInputStream dis) throws IOException
			{
			final CovItem item = new CovItem();
			try {
				item.sample_idx = dis.readInt();
				}
			catch(EOFException err) {
				return null;
				}
			item.pos= dis.readInt();
			item.depth = dis.readFloat();
			item.stddev = dis.readFloat();
			return item;
			}
		@Override
		public void encode(final DataOutputStream dos,final  CovItem o) throws IOException
			{
			dos.writeInt(o.sample_idx);
			dos.writeInt(o.pos);
			dos.writeFloat(o.depth);
			dos.writeFloat(o.stddev);
			}
		@Override
		public CovItemCodec clone()
			{
			return new CovItemCodec();
			}
		}
	

	@Override
	public int doWork(final List<String> args) {
		
		VariantContextWriter w = null;
		try
			{
			final IndexCovUtils indexCovUtils = new IndexCovUtils(this.indexCovTreshold);
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.refPath);
			final SamReaderFactory samReaderFactory = SamReaderFactory.
						makeDefault().
						referenceSequence(CoverageMatrix.this.refPath).
						validationStringency(ValidationStringency.LENIENT)
						;
			
			 final List<Path> inputBams =  IOUtils.unrollPaths(args);
			
			if(inputBams.size()<3) {
				LOG.error("not enough input bam file defined.");
				return -1;
				}
			
			final Set<String> sampleSet = new TreeSet<>();
			final List<String> idx2samples= new ArrayList<String>(inputBams.size());
			for(final Path path: inputBams) {
				try(SamReader sr = samReaderFactory.open(path)) {
					final SAMFileHeader header= sr.getFileHeader();
					
					final String sample = header.getReadGroups().stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtil.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(path));
					if(sampleSet.contains(sample)) {
						LOG.error("duplicate sample "+sample);
						return -1;
						}
					sampleSet.add(sample);
					idx2samples.add(sample);
					}
				}
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			
			w = this.writingVariantsDelegate.dictionary(dict).open(this.outputFile);
			
			final VCFFormatHeaderLine fmtNormDepth = new VCFFormatHeaderLine("D",1,VCFHeaderLineType.Float,"norm Depth");
			metaData.add(fmtNormDepth);
			final VCFFormatHeaderLine fmtStdDev = new VCFFormatHeaderLine("STDDEV",1,VCFHeaderLineType.Float,"standard deviation");
			metaData.add(fmtStdDev);

			final VCFInfoHeaderLine infoStdDev = new VCFInfoHeaderLine(fmtStdDev.getID(),1,VCFHeaderLineType.Float,"standard deviation");
			metaData.add(infoStdDev);
			final VCFInfoHeaderLine infoMedianD = new VCFInfoHeaderLine("MEDIAN",1,VCFHeaderLineType.Float,"median depth");
			metaData.add(infoMedianD);
			final VCFInfoHeaderLine infoNSamples = new VCFInfoHeaderLine("NSAMPLES",1,VCFHeaderLineType.Integer,"number of samples");
			metaData.add(infoNSamples);
			final VCFInfoHeaderLine infoSamples = new VCFInfoHeaderLine("SAMPLES",VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Samples");
			metaData.add(infoSamples);
			final VCFFilterHeaderLine filterAll = new VCFFilterHeaderLine("ALL_AFFECTED","All Samples carry a variant");
			metaData.add(filterAll);
			
			
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true, VCFConstants.END_KEY);
			VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.GENOTYPE_KEY);
			final VCFHeader vcfheader= new VCFHeader(metaData,idx2samples);
			vcfheader.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, vcfheader);
			w.writeHeader(vcfheader);
			
			for(final SAMSequenceRecord ssr: dict.getSequences()) {
				if(!StringUtils.isBlank(restrictContig) && !restrictContig.equals(ssr.getSequenceName())) continue;
				final int depth[]= new int[ssr.getSequenceLength()];
				final BitSet blackListedPositions = new BitSet(depth.length);
								
				// fill black listed regions
				if(this.blackListedPath!=null) {
					try(TabixReader tbr= new TabixReader(this.blackListedPath.toString())) {
						final ContigNameConverter cvt = ContigNameConverter.fromContigSet(tbr.getChromosomes());
						final String ctg = cvt.apply(ssr.getSequenceName());
						if(!StringUtils.isBlank(ctg)) {
							final BedLineCodec codec = new BedLineCodec();
							final TabixReader.Iterator tbxr = tbr.query(ctg,1, ssr.getSequenceLength());
							for(;;) {
								final String line = tbxr.next();
								if(line==null) break;
								final BedLine bed = codec.decode(line);
								if(bed==null) continue;
								int p1 = Math.max(bed.getStart(),1);
								while(p1 <= ssr.getSequenceLength()  && p1 <= bed.getEnd()) {
									blackListedPositions.set(p1-1);
									++p1;
									}
								}
							}
						}
					catch(Throwable err) {
						LOG.warn(err);
						}
					}
				
				final SortingCollection<CovItem> sorter = SortingCollection.newInstance(
						CovItem.class, new CovItemCodec(),
						(A,B)->A.compare0(B),
						this.writingSortingCollection.getMaxRecordsInRam(),
						this.writingSortingCollection.getTmpPaths()
						);
				
				
				for(int bam_idx=0;bam_idx<inputBams.size();++bam_idx) {
					final Path path = inputBams.get(bam_idx);
					LOG.info(ssr.getContig()+":"+path+" "+bam_idx+"/"+inputBams.size());
					try(SamReader sr = samReaderFactory.open(path)) {
						final SAMFileHeader header= sr.getFileHeader();
						
						SequenceUtil.assertSequenceDictionariesEqual(dict,header.getSequenceDictionary());
						Arrays.fill(depth, 0);
						try(CloseableIterator<SAMRecord> siter = sr.queryOverlapping(ssr.getContig(), 1, ssr.getLengthOnReference())) {
							while(siter.hasNext()) {
								final SAMRecord rec= siter.next();
								if(rec.getReadUnmappedFlag()) continue;
								if(!SAMRecordDefaultFilter.accept(rec, this.min_mapq)) continue;
								int ref=rec.getStart();
								final Cigar cigar = rec.getCigar();
								if(cigar==null) continue;
								
								for(CigarElement ce:cigar) {
									final CigarOperator op = ce.getOperator();
									final int len = ce.getLength();
									if(op.consumesReferenceBases()) {
										if(op.consumesReadBases()) {
											for(int i=0;i< len;i++) {
												final int pos = ref+i;
												if(pos < 1) continue;
												if(pos > ssr.getLengthOnReference()) break;
												depth[pos-1]++;
												}
											}
										ref+=len;
										}
									}// loop cigar
								}// end samItere
						} // try
		
					final DiscreteMedian<Integer> discreteMedian = new DiscreteMedian<>();
					int pos=0;
					while(pos< depth.length) {
						if(!blackListedPositions.get(pos) && depth[pos]<=this.max_depth) {
							discreteMedian.add(depth[pos]);
							}
						++pos;
						}
					final double median = discreteMedian.getMedian().orElse(1.0);
					LOG.info(idx2samples.get(bam_idx)+ " :"+ssr.getSequenceName()+" median depth:"+median);
					
					final DiscreteMedian<Integer> localMedian = new DiscreteMedian<>();
					pos=0;
					while(pos< depth.length) {
						if(blackListedPositions.get(pos) /* non pas maxdepth */) {
							++pos;
							continue;
							}
						int pos2=pos;
						localMedian.clear();
						while(pos2 -pos < this.bin_size && pos2< depth.length && !blackListedPositions.get(pos2)) {
							// consider this.max_depth here ?
							localMedian.add(depth[pos2]);
							++pos2;
							}
						if(pos2 -pos == this.bin_size) {
							final double localMed = localMedian.getMedian().orElse(0.0);
							final CovItem item = new CovItem();
							item.pos = pos;
							item.sample_idx = bam_idx;
							item.depth = (float)(localMed/median);
							item.stddev = (float)localMedian.getStandardDeviation().orElse(-1.0);
							sorter.add(item);
							}
						pos = pos2;
						}
					}//end loop over samples
				}//end loop over bams
			sorter.doneAdding();
			sorter.setDestructiveIteration(true);
				
			final CloseableIterator<CovItem> iter = sorter.iterator();
			final EqualRangeIterator<CovItem> iter2 = new EqualRangeIterator<>(iter,(A,B)->Integer.compare(A.pos, B.pos));
			final Allele REF = Allele.create("N", true);
			final Allele DEL = Allele.create("<DEL>", false);
			final Allele DUP = Allele.create("<DUP>", false);
			while(iter2.hasNext()) {
				final List<CovItem> list = iter2.next();
				final CovItem first = list.get(0);
				final double avg_depth = list.stream().mapToDouble(F->F.depth).average().orElse(0);
				final double sum =  list.stream().mapToDouble(F->F.depth).map(D->Math.pow(D - avg_depth,2.0)).sum();
				final double stdDev = Math.sqrt(sum/list.size());
				
				final OptionalDouble optMedianOfmedian = Percentile.median().evaluate(list.stream().mapToDouble(I->I.depth));
				final double medianOfmedian = optMedianOfmedian.orElse(1.0);
				if(medianOfmedian<=0) continue;
				for(int i=0;i< list.size();i++) {
					list.get(i).depth/=medianOfmedian;
					}
				if( list.stream().allMatch(F->Float.isNaN(F.depth) || Float.isInfinite(F.depth))) continue;
				
				final VariantContextBuilder vcb = new VariantContextBuilder();
				vcb.chr(ssr.getContig());
				vcb.start(first.pos+1);
				vcb.stop(first.pos+this.bin_size);
				vcb.attribute(VCFConstants.END_KEY, first.pos+this.bin_size);
				vcb.attribute(infoStdDev.getID(), stdDev);
				vcb.attribute(infoMedianD.getID(), medianOfmedian);
				
				final Set<Allele> alleles = new HashSet<>();
				alleles.add(REF);
				final List<Genotype> genotypes = new ArrayList<>(list.size());
				final Set<String> affected= new TreeSet<>();
			
				for(int i=0;i< list.size();i++) {
					final CovItem item = list.get(i);
					final String sn = idx2samples.get(item.sample_idx);
					final GenotypeBuilder gb;
					switch(indexCovUtils.getType(item.depth))
						{
						case AMBIGOUS: gb = new GenotypeBuilder(sn,Arrays.asList(Allele.NO_CALL,Allele.NO_CALL));break;
						case HET_DEL: alleles.add(DEL);gb = new GenotypeBuilder(sn,Arrays.asList(REF,DEL));affected.add(sn);break;
						case HOM_DEL: alleles.add(DEL);gb = new GenotypeBuilder(sn,Arrays.asList(DEL,DEL));affected.add(sn);break;
						case HET_DUP: alleles.add(DUP);gb = new GenotypeBuilder(sn,Arrays.asList(REF,DUP));affected.add(sn);break;
						case HOM_DUP: alleles.add(DUP);gb = new GenotypeBuilder(sn,Arrays.asList(DUP,DUP));affected.add(sn);break;
						case REF: gb = new GenotypeBuilder(sn,Arrays.asList(REF,REF));break;
						default: throw new IllegalStateException();
						}
					gb.attribute(fmtNormDepth.getID(), item.depth);
					gb.attribute(fmtStdDev.getID(), item.stddev);
					genotypes.add(gb.make());
					}
				if(affected.isEmpty()) continue;
				
				if(affected.size()==inputBams.size()) {
					vcb.filter(filterAll.getID());
					}
				else
					{
					vcb.passFilters();
					}
				vcb.attribute(infoSamples.getID(), new ArrayList<>(affected));
				vcb.attribute(infoNSamples.getID(), affected.size());
				
				vcb.genotypes(genotypes);
				vcb.alleles(alleles);
				w.add(vcb.make());
				}
			iter2.close();
			iter.close();
			sorter.cleanup();
			System.gc();
			}// end while iter
			w.close();w=null;
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			}
		}

public static void main(final String[] args) {
	new CoverageMatrix().instanceMainWithExit(args);
	}

}
