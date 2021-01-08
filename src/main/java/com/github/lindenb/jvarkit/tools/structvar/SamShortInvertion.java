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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Consumer;
import java.util.function.ToIntBiFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.FilteringSamIterator;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.filter.SecondaryOrSupplementaryFilter;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC

## Motivation

finds the regions having some short inversions.

input is a set of BAM files. One file ending with '.list' is interpreted as a file containing some path to the bams.

output is a VCF file

## Example:

```
$ find DIR -type f -name "*.bam" > bams.list
$ java -jar ${JVARKIT_DIST}/samshortinvert.jar -R ref.fasta bams.list |\
 	bcftools view -i 'INFO/DPMAX>10' > out.vcf
```

## Screenshot

* https://twitter.com/yokofakun/status/1222848286048112641
![https://twitter.com/yokofakun/status/1222848286048112641](https://pbs.twimg.com/media/EPhuCJnX4AA3Brc?format=png&name=medium)

* https://twitter.com/yokofakun/status/1222832425518141442
![https://twitter.com/yokofakun/status/1222832425518141442](https://pbs.twimg.com/media/EPhfm8EW4AAiaBq?format=png&name=medium)

* https://twitter.com/yokofakun/status/1222853635656364032
![https://twitter.com/yokofakun/status/1222853635656364032](https://pbs.twimg.com/media/EPhy5fAUYAAMunf?format=png&name=medium)

END_DOC

**/

@Program(name="samshortinvert",
	description="Scan short inversions in SAM using supplementary reads.",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20201210",
	creationDate="20140228"
	)
public class SamShortInvertion extends Launcher
	{
	private static final Logger LOG = Logger.build(SamShortInvertion.class).make();
	private static final byte SUPPORTING_LEFT=(byte)1;
	private static final byte SUPPORTING_RIGHT=(byte)2;
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFaidx = null;
	@Parameter(names={"-m","--maxsize"},description="max size of inversion.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int max_size_inversion = 10_000 ;
	@Parameter(names={"-B","--bed","-r","--rgn"},description=IntervalListProvider.OPT_DESC,splitter= NoSplitter.class,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervallistProvider = null;
	@Parameter(names={"-partition","--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names={"-F","--ratio"},description="Two intervals are the same if they both have more or equals of this fraction of length in common. " + FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private double region_are_same_ratio=0.75;
	@Parameter(names={"-s","-supporting"},description="Don't print the variant if INFO/DP <= 's'")
	private int min_supporting_reads = 1;
	@Parameter(names={"--debug"},description="Debug",hidden=true)
	private boolean debug = false;
	@Parameter(names={"--keep"},description="keep interval if they were used.",hidden=true)
	private boolean keep_flag = false;
	@Parameter(names={"--mapq"},description="min mapping quality")
	private int mapq=1;

	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();

	private static class Arc
		{
		String sample;
		int tid;
		int chromStart;
		int chromEnd;
		boolean consummed=false;
		byte type;
		
		int length() {
			return chromEnd-chromStart+1;
		}
		
		@Override
		public String toString() {
			return "("+tid+"):"+chromStart+"-"+chromEnd +" type="+(int)type;
			}
		}
	
	private boolean testOverlapping2(final int start1,final int end1,final int start2,final int end2 )
		{
		double lenA = CoordMath.getLength(start1, end1);
		double len2 = CoordMath.getOverlap(start1, end1, start2, end2);
		return len2/lenA >= this.region_are_same_ratio;
		}
	
	private boolean testOverlapping(final Locatable a,final Arc b )
		{
		if(!CoordMath.overlaps(a.getStart(), a.getEnd(), b.chromStart, b.chromEnd)) return false;
		return testOverlapping2(a.getStart(), a.getEnd(), b.chromStart, b.chromEnd) && 
				testOverlapping2(b.chromStart, b.chromEnd,a.getStart(), a.getEnd());
		}

	private void dump(
			final SAMSequenceDictionary dict,
			final IntervalTreeMap<List<Arc>> database,
			final VariantContextWriter vcw,
			final Set<String> samples,
			final Integer before
			) {
		if(this.debug) LOG.debug("dump");
		final Allele REF = Allele.create("N", true);
		final Allele SPLIT = Allele.create("<INV>", false);
		final Comparator<Locatable> cmp = new ContigDictComparator(dict).createLocatableComparator();
		final List<SimpleInterval> intervals  = database.keySet().stream().
				filter(R->(before==null?true:R.getEnd() > before.intValue())).
				map(R-> new SimpleInterval(R)).
				sorted(cmp).
				collect(Collectors.toList());
		
		for(final SimpleInterval interval0:intervals) {
			
			final List<Arc> arcs = database.getOverlapping(interval0).
					stream().
					flatMap(L->L.stream()).
					filter(A->!A.consummed).
					filter(A->testOverlapping(interval0, A)).
					collect(Collectors.toList());
			
			if(arcs.isEmpty()) continue;
			if(!keep_flag) arcs.forEach(A->A.consummed=true);
			
			int maxdp = 0;
			final VariantContextBuilder vcb = new VariantContextBuilder();
			final Set<Allele> alleles = new HashSet<>();
			alleles.add(REF);
			final List<Genotype> genotypes = new ArrayList<>(samples.size());
			
			
			vcb.chr(dict.getSequence(arcs.get(0).tid).getSequenceName());
			final int chromStart = arcs.stream().mapToInt(A->A.chromStart).min().getAsInt();
			vcb.start(chromStart);
			final int chromEnd = arcs.stream().mapToInt(A->A.chromEnd).max().getAsInt();
			vcb.stop(chromEnd);
			
			vcb.attribute(VCFConstants.END_KEY, chromEnd);
			vcb.attribute("SVLEN", CoordMath.getLength(chromStart,chromEnd));
			
			int depth = 0;
			final Set<String> gotSamples = new TreeSet<>();
			for(final String sample : samples) {
			
				final List<Arc> sampleArcs = arcs.stream().
						filter(A->A.sample.equals(sample)).
						collect(Collectors.toList());
				if(sampleArcs.isEmpty())
					{
					genotypes.add(GenotypeBuilder.createMissing(sample, 2));
					}
				else
					{
					final GenotypeBuilder gb = new GenotypeBuilder(sample);
					alleles.add(SPLIT);
					gb.alleles(Arrays.asList(REF,SPLIT));
					final int countCat1= (int)sampleArcs.stream().filter(A->A.type==SUPPORTING_LEFT).count();
					final int countCat2= (int)sampleArcs.stream().filter(A->A.type==SUPPORTING_RIGHT).count();
					gb.DP(countCat1+countCat2);
					gb.attribute("N5", countCat1);
					gb.attribute("N3", countCat2);
					
					maxdp = Math.max(maxdp, countCat1+countCat2);
					
					depth+=countCat1+countCat2;
					genotypes.add(gb.make());
					gotSamples.add(sample);
					}
				}
			if(depth<=this.min_supporting_reads) continue;
			
			vcb.genotypes(genotypes);
			vcb.alleles(alleles);
			vcb.attribute(VCFConstants.DEPTH_KEY, depth);
			vcb.attribute("NSAMPLES", gotSamples.size());
			vcb.attribute("SAMPLES", new ArrayList<>(gotSamples));
			vcb.attribute(VCFConstants.SVTYPE, "INV");
			vcb.attribute("DPMAX", maxdp);
			vcw.add(vcb.make());
			}
		
		}

	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_size_inversion<100) {
			LOG.error("max size insersion must be >=100");
			return -1;
			}
	 	final Map<SamReader,CloseableIterator<SAMRecord>> samReaders = new HashMap<>();
	 	VariantContextWriter vcw= null;
	 	final  IntervalTreeMap<List<Arc>> database = new IntervalTreeMap<>(); 
		try {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceFaidx);
			
			final short SA_TAG = SAMTag.SA.getBinaryTag();

			
			final QueryInterval queryIntervals[] = this.intervallistProvider==null?
				null:
				this.intervallistProvider.
					dictionary(dict).
					optimizedQueryIntervals()
					;
			
			final AggregateFilter theFilter = new AggregateFilter(
					Arrays.asList(
							new MappingQualityFilter(this.mapq),
							new DuplicateReadFilter(),
							new SecondaryOrSupplementaryFilter(),
							new FailsVendorReadQualityFilter(),
							new SamRecordFilter() {
								@Override
								public boolean filterOut(SAMRecord first, SAMRecord second) {
									return filterOut(first) || filterOut(second);
								}
								@Override
								public boolean filterOut(final SAMRecord rec) {
									if(rec.getReadUnmappedFlag()) return true;
									if(rec.getAttribute(SA_TAG) == null) return true;
									final Cigar cigar = rec.getCigar();
									if(cigar==null || cigar.isEmpty() || !cigar.isClipped()) return true;
									return false;
								}
							}
							)				
					);
			
			for(final Path samPath:IOUtils.unrollPaths(args)) {
				final SamReader srf = SamReaderFactory.
						makeDefault().
						validationStringency(ValidationStringency.LENIENT).
						referenceSequence(this.referenceFaidx).
						open(samPath);
				
				final CloseableIterator<SAMRecord> iter;
				if(queryIntervals!=null)
					{
					iter = srf.query(queryIntervals,false);
					}
				else
					{
					iter = srf.iterator();
					}
				final FilteringSamIterator sfi = new FilteringSamIterator(iter,theFilter);
				samReaders.put(srf,sfi);
				}
			final SamFileHeaderMerger headerMerger  = new SamFileHeaderMerger(SAMFileHeader.SortOrder.coordinate,samReaders.keySet().stream().map(SR->SR.getFileHeader()).collect(Collectors.toList()),false);
			final MergingSamRecordIterator iter = new MergingSamRecordIterator( headerMerger,samReaders, true);
			
			
			
			final Set<String> samples = headerMerger.getHeaders().stream().
					flatMap(R->R.getReadGroups().stream()).
					map(RG->this.partition.apply(RG, null)).
					filter(S->!StringUtil.isBlank(S)).
					collect(Collectors.toCollection(TreeSet::new));
			
			if(samples.isEmpty())
				{
				iter.close();
				LOG.error("No samples/bam defined");
				return -1;
				}
			
			final ToIntBiFunction<Locatable, Locatable> distance = (A,B) ->{
				if(CoordMath.overlaps(A.getStart(), A.getEnd(), B.getStart(), B.getEnd())) return 0;
				if(A.getEnd()<B.getStart()) {
					return B.getStart() - A.getEnd();
					}
				else
					{
					return A.getStart() - B.getEnd();
					}
				};
		
			
			final Set<VCFHeaderLine> meta=new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(meta,true,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.DEPTH_KEY
					);
			VCFStandardHeaderLines.addStandardInfoLines(meta,true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					
					);
			meta.add(new VCFFormatHeaderLine("N5", 1, VCFHeaderLineType.Integer,"Number of validating clipped reads in 5'"));
			meta.add(new VCFFormatHeaderLine("N3", 1, VCFHeaderLineType.Integer,"Number of validating clipped reads in 3'"));
			meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
			meta.add(new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of sample having some split reads"));
			meta.add(new VCFInfoHeaderLine("SAMPLES", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String,"Samples having some split reads"));
			meta.add(new VCFInfoHeaderLine("DPMAX", 1, VCFHeaderLineType.Integer,"MAX DP among samples"));
			meta.add(new VCFInfoHeaderLine("SVTYPE", 1, VCFHeaderLineType.String,"Structural variant type"));
			
			
			
			final VCFHeader header=new VCFHeader(meta,samples);
			JVarkitVersion.getInstance().addMetaData(this, header);
			header.setSequenceDictionary(dict);
			vcw = this.writingVariants.open(this.outputFile);
			vcw.writeHeader(header);
			
			
			final ProgressFactory.Watcher<SAMRecord> progress= ProgressFactory.
					newInstance().
					dictionary(dict).
					logger(LOG).
					build();
			
			String prevContig=null;
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.apply(iter.next());
				
				if(theFilter.filterOut(rec)) continue;
				
				final String sample= this.partition.getPartion(rec, null);
				if(StringUtil.isBlank(sample))continue;
				
				final List<SAMRecord> others = SAMUtils.getOtherCanonicalAlignments(rec).
						stream().
						filter(R->rec.getContig().equals(R.getContig())).
						filter(R->rec.getReadNegativeStrandFlag()!=R.getReadNegativeStrandFlag()).
						filter(R->distance.applyAsInt(rec,R)< this.max_size_inversion).
						collect(Collectors.toList());
				
				if(others.isEmpty()) continue;
				
				if(!rec.getContig().equals(prevContig)) {
					dump(dict,database,vcw,samples,null);
					database.clear();
					prevContig = rec.getContig();
					}
				else
					{
					final int before = rec.getUnclippedStart() - this.max_size_inversion*2;
					dump(dict,database,vcw,samples,before);
					database.entrySet().removeIf(entries->entries.getKey().getEnd()< before);
					}

				final Consumer<Arc> registerArc = (A)->{
					if(A.chromEnd<=A.chromStart) throw new IllegalArgumentException(A.toString());
					final Interval rgn = new Interval(rec.getContig(), A.chromStart,A.chromEnd);
					List<Arc> list = database.get(rgn);
					if(list==null) {
						list = new ArrayList<>();
						database.put(rgn,list);
						}
					list.add(A);					
					};
				final Cigar cigar = rec.getCigar();
				if(cigar.isLeftClipped())
					{
					for(final SAMRecord rec2:others) {
						// NON if(rec.getEnd()>= rec2.getStart()) continue;
						final Arc arc = new Arc();
						arc.sample  = sample;
						arc.tid = rec.getReferenceIndex();
						arc.chromStart = Math.min(rec.getStart(),rec2.getStart());
						arc.chromEnd = Math.max(rec.getEnd(),rec2.getEnd());
						if(arc.length()> this.max_size_inversion) continue;

						
						arc.type = SUPPORTING_LEFT;
						registerArc.accept(arc);
						}
					}
			
				
				if(cigar.isRightClipped())
					{
					for(final SAMRecord rec2:others) {
						final Arc arc = new Arc();
						arc.sample  = sample;
						arc.tid = rec.getReferenceIndex();
						arc.chromStart = Math.min(rec.getStart(),rec2.getStart());
						arc.chromEnd = Math.max(rec.getEnd(),rec2.getEnd());
						if(arc.length()> this.max_size_inversion) continue;
						
						arc.type = SUPPORTING_RIGHT;
						registerArc.accept(arc);
						}
					}				
				}
			dump(dict,database,vcw,samples,null);
			iter.close();
			progress.close();
			vcw.close();vcw=null;
			for(final SamReader sr: samReaders.keySet()) {
				samReaders.get(sr).close();
				sr.close();
			}
			return 0;
			}
		catch (final Throwable e) {
			LOG.error(e);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcw);	
			}
		}

	
	public static void main(final String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
