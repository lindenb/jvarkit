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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.Decoy;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
/**

BEGIN_DOC


```
$ java -jar dist/svneg.jar --bams jeter.list src/test/resources/HG02260.transloc.chr9.14.bam

##fileformat=VCFv4.2
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=SR,Number=4,Type=Integer,Description="Supporting reads: contig1-forward,contig1-reverse,contig2-forward,contig2-reverse">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=CHROM2,Number=1,Type=String,Description="other chromosome">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##INFO=<ID=POS2,Number=1,Type=Integer,Description="other position">
##INFO=<ID=STDDEV_POS1,Number=1,Type=Integer,Description="std deviation to position 1">
##INFO=<ID=STDDEV_POS2,Number=1,Type=Integer,Description="std deviation to position 2">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variation type">
##svneg.meta=compilation:20190415102101 githash:4c6799f7 htsjdk:2.19.0 date:20190415102110 cmd:--bams jeter.list src/test/resources/HG02260.transloc.chr9.14.bam
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG02260
9	137230996	9:137230996:14:79839048	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=14;DP=18;POS2=79839048;STDDEV_POS1=120;STDDEV_POS2=187;SVTYPE=BND	GT:DP:SR	0/1:18:5,13,13,5
14	79839119	14:79839119:9:137230968	N	<TRANSLOC>	18	.	AC=1;AF=0.500;AN=2;CHROM2=9;DP=18;POS2=137230968;STDDEV_POS1=154;STDDEV_POS2=145;SVTYPE=BND	GT:DP:SR	0/1:18:13,5,5,13


```

END_DOC
*/
@Program(name="svneg",
	description="Find Structural Variation by Negative Comparaison",
	keywords={"sam","bam","sv","translocation"},
	creationDate="20190413",
	modificationDate="20190826"
	)
public class SamStructVarNeg extends Launcher {
	private static final Logger LOG = Logger.build(SamStructVarNeg.class).make();
	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description="For Reading CRAM. " + INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path refPath = null;
	@Parameter(names={"-c","--min-cigar-size"},description="Min cigar size for clipped cigar elements.")
	private int min_cigar_length = 10;
	@Parameter(names={"-b","--controls"},description="Control Bams. a fil ending with '.list' is interpreted as a list of path, one per line.",required=true)
	private List<Path> controlBamPaths = new ArrayList<>();
	@Parameter(names={"--mapq"},description="min mapping quality.")
	private int min_mapq = 0;
	@Parameter(names={"-x","--exclude"},description="Optional BED file Excluding. SV shouldn't overlap this bed.")
	private Path excludeBedFile = null;
	
	private SamReaderFactory samReaderFactory = null;
	private final Decoy decoy = Decoy.getDefaultInstance();
	
	
	private abstract interface SVChecker {
		boolean check(final SAMRecord sr);
	}
	
	private abstract class TranslocChecker implements SVChecker {
		String ctg1;
		String ctg2;
		
	}
	
	private static class Window extends SimpleInterval{
		Window(final Locatable loc) {
			super(loc);
		}
		
	}
	
	private static interface SVScanner {
		Locatable getInterval();
		Predicate<SAMRecord> getPredicate();
		VariantContext makeVariant();
		}
	
	private class SplitLoc {
		final String contig;
		final int pos;
		SplitLoc(final String contig,int pos) {
			this.contig = contig;
			this.pos = pos;
			}
		boolean withinDistanceOf(final SplitLoc o,final int d) {
			return contig.equals(o.contig) && 
					Math.abs(this.pos-o.pos)< d;
			}
		}
	
	private class SplitPair {
		private final SplitLoc split1;
		private final SplitLoc split2;
		SplitPair( final SplitLoc loc1) {
			this(loc1,loc1);
			}
		SplitPair( final SplitLoc loc1, final SplitLoc loc2) {
			this.split1 = loc1;
			this.split2 = loc2;
			}
	
		}
	
	private void getPredicates(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return Collections.emptyList();
		if(rec.getReadPairedFlag() && 
			!rec.getMateUnmappedFlag()) 
			{
			final String ctg = rec.getMateReferenceName();
			
			if(ctg.equals(rec.getContig()) && 
				!rec.getProperPairFlag() && 
				rec.getEnd() < rec.getMateAlignmentStart()
				) {
				final Predicate<SAMRecord> pred=REC->{
					if(REC.getReadUnmappedFlag()) return false;
					if(!REC.getReadPairedFlag()) return false;
					if(REC.getProperPairFlag()) return false;
					if(!REC.getContig().equals(rec.getContig())) return false;
					if(!REC.getMateReferenceName().equals(rec.getMateReferenceName())) return false;
					if(!(REC.getEnd()<REC.getMateAlignmentStart())) return false;
					if(!withinDistance(rec.getEnd(),REC.getEnd()),50) return false;
					if(!withinDistance(rec.getMateAlignmentStart(),REC.getMateAlignmentStart()),50) return false;
					return true;
					};
				}
			else if(!this.decoy.isDecoy(ctg))
				{
				
				}
			}
		final Predicate<SAMRecord> cigarPredicate = (SR)->{
			final Cigar cigar = SR.getCigar();
			if(cigar==null || cigar.isEmpty()) return false;
			if(cigar.numCigarElements()<=1) return false;
			for(int side=0;side<2;side++) {
				final CigarElement ce = side==0?cigar.getFirstCigarElement():cigar.getLastCigarElement();
				if(!ce.getOperator().isClipping()) continue;
				if(ce.getLength()> min_cigar_length) return true;
				}
			return false;
			};
		
		if(!cigarPredicate.test(rec)) return false;	
			
		for(final SAMRecord rec2:SAMUtils.getOtherCanonicalAlignments(rec)) {
			final String ctg = rec2.getContig();
			if(this.decoy.isDecoy(ctg)) continue;
			
			final Predicate<SAMRecord> pred=REC->{
				if(REC.getReadUnmappedFlag()) return false;
				if(!REC.getContig().equals(rec.getContig())) return false;
				if(!cigarPredicate.test(REC)) return false;
				
				final Cigar cigar2 = REC.getCigar();
				if(cigar2.numCigarElements()<=1) return false;

				
				for(final SAMRecord REC2:SAMUtils.getOtherCanonicalAlignments(REC))
					{
					if(!REC2.getContig().equals(rec2.getContig())) return false;
					
					
					if(!withinDistance(rec.getMateAlignmentStart(),REC.getMateAlignmentStart()),50) return false;
					
					}
				return false;
				};
			}
		}
	
	private boolean accept(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return false;
		if(rec.getMappingQuality() < this.min_mapq) return false;
		if(rec.getDuplicateReadFlag()) return false;
		if(rec.isSecondaryOrSupplementary()) return false;
		if(rec.getReadFailsVendorQualityCheckFlag()) return false;
		if(this.decoy.isDecoy(rec.getContig())) return false;
		if(rec.getReadPairedFlag()) {
			if(!rec.getMateUnmappedFlag()) {
				if(this.decoy.isDecoy(rec.getMateReferenceName())) return false;
				}
			}
		return true;
		}
	
	private SVScanner getSVScanner(final SAMRecord rec) {
		if(!accept(rec)) return null;
		final Cigar cigar = rec.getCigar();
		if(cigar.numCigarElements()>1) {
			for(int side=0;side<2;++side) {
				final CigarElement ce = (side==0?cigar.getFirstCigarElement():cigar.getLastCigarElement());
				if(!ce.getOperator().isClipping()) continue;
				if(ce.getLength()<this.min_cigar_length) continue;
				final Locatable rgn = (side==0?rec:rec);
				return new SVScanner() {
					public Locatable getInterval() 
						{
						return rgn;
						}
					@Override
					public Predicate<SAMRecord> getPredicate() {
						return null;
						}
					@Override
					public VariantContext makeVariant() {
						return null;
						}
					};
				}
			}
		
		if(rec.getReadPairedFlag()) {
			if(rec.getMateUnmappedFlag()) return null;
			if(rec.getProperPairFlag()) return null;
			if(!rec.getReferenceIndex().equals(rec.getMateReferenceIndex())) {
				return new SVScanner() {
					@Override
					public VariantContext makeVariant() {
						final Allele REF=Allele.create("N", true);
						final Allele ALT=Allele.create("<TRANSLOC>", false);
						final VariantContextBuilder vcb = new VariantContextBuilder(null, rec.getContig(), rec.getStart(),rec.getEnd(), Arrays.asList(REF,ALT));
						return vcb.make();
						}
					@Override
					public Locatable getInterval() {
						return rec;
						}
					@Override
					public Predicate<SAMRecord> getPredicate() {
						return (SR)->SR.getContig().equals(rec.getContig()) &&
								SR.getReadPairedFlag() && 
								!SR.getMateUnmappedFlag() &&
								SR.getMateReferenceName().equals(rec.getMateReferenceName());
						}
					};
				}
			}
		return null;
		}
	
	private void recursive(int case_idx,final List<SamReader> casesReaders,final Window win,final List<Path> controlBams,final VariantContextWriter vcw) throws IOException {
		final Predicate<SAMRecord> predicate = param.getPredicate();
		if(case_idx==casesReaders.size()) {
			int count_controls = 0;
			for(final Path ctrlPath: controlBams) {
				try(SamReader sr = this.samReaderFactory.open(ctrlPath)) {
					try(CloseableIterator<SAMRecord> iter = sr.queryOverlapping(rgn.getContig(), rgn.getStart(),rgn.getEnd())) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							if(!accept(rec)) continue;
							if(predicate.test(iter.next())) {
								count_controls++;
								break;
								}
							}
						}
					}
				if(count_controls>0) return;
				}
			// TODO run report !
			return;
			}
		boolean got_sv = false;
		try(CloseableIterator<SAMRecord> iter = casesReaders.get(case_idx).queryOverlapping(win.getContig(), win.getStart(),win.getEnd())) {
			while(iter.hasNext()) {
				if(predicate.test(iter.next())) {
					got_sv=true;
					break;
				}
			}
		}
		if(!got_sv) return;
		recursive(case_idx+1,casesReaders,win,controlBams,vcw);
		}
	private void extractSlipt(final SAMRecord rec1) {
		final Cigar cigar1 = rec1.getCigar();
		if(cigar1==null || cigar1.isEmpty() || cigar1.numCigarElements()<=1) {
			return false;
			}
		final List<SAMRecord> saList = SAMUtils.getOtherCanonicalAlignments(rec1);

		if(rec1.getReadPairedFlag() && 
			!rec1.getMateUnmappedFlag() &&
			!this.decoy.isDecoy(rec1.getContig()) &&
			!rec1.getContig().equals(rec1.getMateReferenceName())) {

		
			}
		
		final List<SAMRecord> saList = SAMUtils.getOtherCanonicalAlignments(rec1);
		
		for(int side=0;side<2;side++) {
			final CigarElement ce = side==0?cigar1.getFirstCigarElement():cigar1.getLastCigarElement();
			if(!ce.getOperator().isClipping()) continue;
			if(ce.getLength() < min_cigar_length) continue;
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		VariantContextWriter out = null;
		
		final List<Path> casePaths = IOUtils.unrollPaths(args);
		
		if(casePaths.isEmpty()) {
			LOG.error("controlPaths illegal number of arguments");
			return -1;
			}
		final List<SamReader> casesReaders = new ArrayList<>(casePaths.size());

		try {
			SAMSequenceDictionary dict = null;
			
			final List<Path> controlPaths;
			if(this.controlBamPaths.size()==1 && this.controlBamPaths.get(0).toString().endsWith(".list")) {
				controlPaths = IOUtils.unrollPath(this.controlBamPaths.get(0));
			} else {
				controlPaths = new ArrayList<>(this.controlBamPaths);
			}
			
			this.samReaderFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
			if(this.refPath!=null) {
				this.samReaderFactory.referenceSequence(this.refPath);
				dict = SequenceDictionaryUtils.extractRequired(this.refPath);
				}
			
			final IntervalTreeMap<Interval> excludeBedMap ;

			if(this.excludeBedFile!=null) {
				excludeBedMap = new IntervalTreeMap<Interval>();
				final BedLineCodec bedCodec = new BedLineCodec();
				try(BufferedReader br=com.github.lindenb.jvarkit.io.IOUtils.openPathForBufferedReading(this.excludeBedFile)) {
					br.lines().
					filter(line->!(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty())).
					map(line->bedCodec.decode(line)).
					filter(B->B!=null).
					map(B->B.toInterval()).
					filter(L->L.getStart()<L.getEnd()).
					forEach(B->{
						excludeBedMap.put(B,B);							
						});	
					}
				}
			else
				{
				excludeBedMap = null;
				}
			
			for(final Path p: casePaths) {
				final SamReader sr = this.samReaderFactory.open(p);
				if(!sr.hasIndex())  {
					LOG.error("BAM file is not indexed : "+p);
					sr.close();
					return -1;
					}
				final SAMSequenceDictionary dict2 = SequenceDictionaryUtils.extractRequired(sr.getFileHeader());
				if(dict==null) {
					dict=dict2;
					}
				else if(!SequenceUtil.areSequenceDictionariesEqual(dict, dict2))
					{
					LOG.error(JvarkitException.DictionariesAreNotTheSame.getMessage(dict, dict2));
					return -1;
					}
				casesReaders.add(sr);
				}
			
			final Set<VCFHeaderLine> metaData=new HashSet<>();
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
			final VCFFormatHeaderLine supportingReadsFormat = new VCFFormatHeaderLine("SR",
					4,
					VCFHeaderLineType.Integer,
					"Supporting reads: contig1-forward,contig1-reverse,contig2-forward,contig2-reverse"
					);
			metaData.add(supportingReadsFormat);
			final VCFInfoHeaderLine stdDevContig1Info = new VCFInfoHeaderLine("STDDEV_POS1",
					1,
					VCFHeaderLineType.Integer,
					"std deviation to position 1"
					);
			metaData.add(stdDevContig1Info);
			final VCFInfoHeaderLine stdDevContig2Info = new VCFInfoHeaderLine("STDDEV_POS2",
					1,
					VCFHeaderLineType.Integer,
					"std deviation to position 2"
					);
			metaData.add(stdDevContig2Info);
			final VCFInfoHeaderLine chrom2Info = new VCFInfoHeaderLine("CHROM2",
					1,
					VCFHeaderLineType.String,
					"other chromosome"
					);
			metaData.add(chrom2Info);
			final VCFInfoHeaderLine pos2Info = new VCFInfoHeaderLine("POS2",
					1,
					VCFHeaderLineType.Integer,
					"other position"
					);
			metaData.add(pos2Info);
			
			final VCFHeader vcfHeader= new VCFHeader(metaData);
			vcfHeader.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, vcfHeader);
			
			out = VCFUtils.createVariantContextWriterToPath(this.outputFile);
			out.writeHeader(vcfHeader);

			
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
					newInstance().
					dictionary(dict).
					logger(LOG).
					build();

			final SamReader first = casesReaders.get(0);
			final List<Window> windows = new ArrayList<>();
			final Map<Interval,Window> winwin = new HashMap<>();
			try(CloseableIterator<SAMRecord> iter = first.iterator()) {
				for(;;) {
					final int window_size=1000;
					final int window_shift=500;
					final SAMRecord rec = iter.hasNext()?progress.apply(iter.next()):null;
					if(rec!=null && !accept(rec)) continue;
					
					int idx=0;
					while(idx < windows.size()) {
						final Window w = windows.get(idx);
						if(w.getContig().equals(rec.getContig()) && 
							w.getEnd() < rec.getUnclippedStart())  {
							winwin.remove(w);
							windows.remove(idx);
							recursive(1,casesReaders,w,controlPaths,out);
							}
						else
							{
							++idx;
							}
						}
					
					if(rec==null || (windows.size()>0 && !windows.get(0).getContig().equals(rec.getContig())))
						{
						for(final Window w:windows) {
							recursive(1,casesReaders,w,controlPaths,out);
							}
						if(rec==null) break;
						windows.clear();
						winwin.clear();
						}
					int x = rec.getStart()-rec.getStart()%window_shift;
					while(x <=rec.getEnd()) {
						final Interval wr =  new Interval(rec.getContig(),x,x+window_size);
						Window w= winwin.get(wr);
						if(w==null) {
							w=new Window(wr);
							winwin.put(wr,w);
							windows.add(w);
							}
						x+=window_shift;
						}
					
					final SVScanner param = this.getSVScanner(rec);
					if(param==null) continue;
					if(excludeBedMap!=null && excludeBedMap.containsOverlapping(rec)) continue;
					}
			}
			
		out.close();
		out=null;
		progress.close();
		
		return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			casesReaders.stream().forEach(S->CloserUtil.close(S));
			CloserUtil.close(out);
			}
		}
	
	public static void main(final String[] args) {
		new SamStructVarNeg().instanceMainWithExit(args);
		}

}
