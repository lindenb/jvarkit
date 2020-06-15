/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
import java.util.Arrays;import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalInt;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.SAMRecordPartition;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.MergingSamRecordIterator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamFileHeaderMerger;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.AggregateFilter;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.FailsVendorReadQualityFilter;
import htsjdk.samtools.filter.FilteringSamIterator;
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
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**
BEGIN_DOC


END_DOC

**/

@Program(name="saminversions",
	description="Scan inversions in SAM",
	keywords={"sam","bam","sv","inversion"},
	modificationDate="20200608",
	creationDate="20200609",
	generate_doc=false
	)
public class SamInversions extends Launcher
	{
	private static final Logger LOG = Logger.build(SamInversions.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-t","--treshold"},description="Two SV are the same if they overlap and share the 'x' fraction of their length. "+ FractionConverter.OPT_DESC,converter=FractionConverter.class)
	private double fraction=0.70;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referenceFaidx = null;
	@Parameter(names={"--max-size"},description="max size of inversion.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int max_size_inversion = 10_000_000 ;
	@Parameter(names={"--min-size"},description="min size of inversion.",splitter=NoSplitter.class,converter=DistanceParser.StringConverter.class)
	private int min_size_inversion = 150 ;
	@Parameter(names={"-B","--bed","-r","--rgn"},description=IntervalListProvider.OPT_DESC,splitter= NoSplitter.class,converter=IntervalListProvider.StringConverter.class)
	private IntervalListProvider intervallistProvider = null;
	@Parameter(names={"-partition","--partition"},description=SAMRecordPartition.OPT_DESC)
	private SAMRecordPartition partition = SAMRecordPartition.sample;
	@Parameter(names={"-s","-supporting"},description="Don't print the variant if max FORMAT/DP < 's'")
	private int min_supporting_reads = 1;
	@Parameter(names={"--mapq"},description="min MAPQ")
	private int mapq = 1;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariants = new WritingVariantsDelegate();

	private static class Arc extends SimpleInterval
		{
		final Counter<String> samples = new Counter<>();
		Arc(final Locatable loc) {
			super(loc);
			}
		}
	
	

	private void dump(
			final SAMSequenceDictionary dict,
			final IntervalTreeMap<Arc> database,
			final VariantContextWriter vcw,
			final Set<String> samples,
			final OptionalInt before
			) {
		final Allele REF = Allele.create("N", true);
		final Allele SPLIT = Allele.create("<INV>", false);
		
		for(Arc arc:database.values()) {
			if(before.isPresent() && arc.getEnd() >= before.getAsInt()) continue;
			final int maxdp = (int)arc.samples.getMaxCount().orElse(0L);
			if(maxdp < this.min_supporting_reads) continue;
			
			
			final VariantContextBuilder vcb = new VariantContextBuilder();
			final Set<Allele> alleles = new HashSet<>();
			alleles.add(REF);
			final List<Genotype> genotypes = new ArrayList<>(samples.size());
			
			
			vcb.chr(arc.getContig());
			final int chromStart = arc.getStart();
			vcb.start(chromStart);
			final int chromEnd =  arc.getEnd();
			vcb.stop(chromEnd);
			
			vcb.attribute(VCFConstants.END_KEY, chromEnd);
			
			int depth = 0;
			int nsamples = 0;
			for(final String sample : samples) {
				if(arc.samples.count(sample)==0L)
					{
					genotypes.add(GenotypeBuilder.createMissing(sample, 2));
					}
				else
					{
					final GenotypeBuilder gb = new GenotypeBuilder(sample);
					alleles.add(SPLIT);
					gb.alleles(Arrays.asList(REF,SPLIT));
					final int dp = (int)arc.samples.count(sample);
					gb.DP(dp);					
					depth+=dp;
					genotypes.add(gb.make());
					++nsamples;
					}
				}
			
			vcb.genotypes(genotypes);
			vcb.alleles(alleles);
			vcb.attribute(VCFConstants.DEPTH_KEY, depth);
			vcb.attribute("NSAMPLES", nsamples);
			vcb.attribute(VCFConstants.SVTYPE, "INV");
			vcb.attribute("SVLEN", CoordMath.getLength(chromStart,chromEnd));
			vcb.attribute("DPMAX", maxdp);
			vcw.add(vcb.make());
			}
		
		}

	
	@Override
	public int doWork(final List<String> args) {
		if(this.max_size_inversion<100) {
			LOG.error("max size inversion must be >=100");
			return -1;
			}
		if(this.fraction<=0.0) {
			LOG.error("bad fraction: " + this.fraction);
			return -1;
			}
	 	final Map<SamReader,CloseableIterator<SAMRecord>> samReaders = new HashMap<>();
	 	VariantContextWriter vcw= null;
	 	final  IntervalTreeMap<Arc> database = new IntervalTreeMap<>(); 
		try {
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(this.referenceFaidx);
			
			final ToIntFunction<SAMRecord> mateEnd = REC->SAMUtils.getMateCigar(REC)!=null?
				SAMUtils.getMateAlignmentEnd(REC):
				REC.getMateAlignmentStart();
			
			final ToIntFunction<SAMRecord> recordLen = REC->{
				final int start = Math.min(REC.getStart(), REC.getMateAlignmentStart());
				final int end = Math.max(REC.getEnd(), mateEnd.applyAsInt(REC));
				return CoordMath.getLength(start, end);
				};
				
			
			final QueryInterval queryIntervals[] = this.intervallistProvider==null?
				null:
				this.intervallistProvider.
					dictionary(dict).
					optimizedQueryIntervals()
					;
			
			final AggregateFilter theFilter = new AggregateFilter(
					Arrays.asList(
							new DuplicateReadFilter(),
							new SecondaryOrSupplementaryFilter(),
							new FailsVendorReadQualityFilter(),
							new SamRecordFilter() {
								@Override
								public boolean filterOut(final SAMRecord first, final SAMRecord second) {
									return filterOut(first) || filterOut(second);
								}
								@Override
								public boolean filterOut(final SAMRecord rec) {
									if(rec.getReadUnmappedFlag()) return true;
									if(rec.getMappingQuality()<mapq) return true;
									if(!rec.getReadPairedFlag()) return true;
 									if(rec.getMateUnmappedFlag()) return true;
									if(!rec.getReferenceIndex().equals(rec.getReferenceIndex())) return true;
									if(rec.getReadNegativeStrandFlag()!=rec.getMateNegativeStrandFlag()) return true;
									final int len =  recordLen.applyAsInt(rec);
									if(len < min_size_inversion) return true;
									if(len > max_size_inversion) return true;
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
			
		
			
			final Set<VCFHeaderLine> meta=new HashSet<>();
			VCFStandardHeaderLines.addStandardFormatLines(meta,true,
					VCFConstants.GENOTYPE_KEY,
					VCFConstants.DEPTH_KEY
					
					);
			VCFStandardHeaderLines.addStandardInfoLines(meta,true,
					VCFConstants.DEPTH_KEY,
					VCFConstants.END_KEY
					);
			meta.add(new VCFInfoHeaderLine("SVLEN", 1, VCFHeaderLineType.Integer,"SV length"));
			meta.add(new VCFInfoHeaderLine("NSAMPLES", 1, VCFHeaderLineType.Integer,"Number of samples with SV."));
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
				
				final SimpleInterval r = new SimpleInterval(rec.getContig() ,
						Math.min(rec.getStart(),rec.getMateAlignmentStart()),
						Math.max(rec.getEnd(),mateEnd.applyAsInt(rec))
						);
				if(r.getLengthOnReference() < min_size_inversion) continue;
				if(r.getLengthOnReference() > max_size_inversion) continue;
				
				if(!rec.getContig().equals(prevContig)) {
					dump(dict,database,vcw,samples, OptionalInt.empty());
					database.clear();
					prevContig = rec.getContig();
					}
				else
					{
					dump(dict,database,vcw,samples,OptionalInt.of(rec.getStart()-2*max_size_inversion));
					database.entrySet().removeIf(KV->KV.getKey().getEnd()+1+(int)Math.ceil(KV.getKey().getLengthOnReference()*(1.0-fraction)) < rec.getStart());
					}
				boolean found=false;
				
				final String sample= this.partition.getPartion(rec,null);
				if(StringUtil.isBlank(sample))continue;

				
				for(Arc arc:database.getOverlapping(r)) {
					final double L1 = arc.getLengthOnReference();
					final double L2 = r.getLengthOnReference();
					final double L3 = r.getIntersectionLength(arc);
					if(L3/L1 < this.fraction ) continue;
					if(L3/L2 < this.fraction ) continue;
					
					arc.samples.incr(sample);
					found=true;
					}
				if(!found) {
					Arc arc = new Arc(r);
					arc.samples.incr(sample);
					database.put(new Interval(arc),arc);
					}
				}
			dump(dict,database,vcw,samples,OptionalInt.empty());
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
		new SamInversions().instanceMainWithExit(args);
		}
	}
