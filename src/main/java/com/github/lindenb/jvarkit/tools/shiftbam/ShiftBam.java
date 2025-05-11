package com.github.lindenb.jvarkit.tools.shiftbam;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.samtools.util.SimplePosition;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC



shift coordinates of a BAM.

## Usage 1

Say the bam has been mapped on a sub-fasta.

```
# here all contigs look like "chrxxx:123-456"
samtools faidx ref.fa "chr1:2345-6789" >  ref2.fa
samtools faidx ref.fa "chr21:12345-16789" >>  ref2.fa




(...)
bwa mem ref2.fa read1.fq read2.fq |\
	java -jar dist/jvarkit.jar shiftbam -R2 ref.Fa

````


In this sub-fasta all the chromosomes NAMES **MUST** look like 'chr1:123-456'.
This tool takes as input the original REF and the bam and it's coordinates in the new ref
and shit the coordinates to the new reference, including the SA:Z attributes.
Ouput BAM is not sorted on coordinate.


## Usage 2

we're creating a mini-genome containing the gene of interest using ROI.bed

```
java -jar dist/jvarkit.jar shiftbam --bed ROI.bed in.bam
```

## Example:

the following Makefile extract a BAM to a sub-region, call delly on a smaller REF, convert back the vcf to the original ref

```makefile
OUTDIR=work
PREFIX=myexperiment.
REF=ref.fa
BAM=${HOME}/jeter.bam
INTERVAL=chr1:123-456

$(OUTDIR)/$(PREFIX).shift.vcf: $(OUTDIR)/$(PREFIX).call.bcf
	bcftools view $< | java -jar $${JVARKIT_DIST}/jvarkit.jar shiftvcf -R2 $(REF) > $@

$(OUTDIR)/$(PREFIX).call.bcf : $(OUTDIR)/shift.01.bam $(OUTDIR)/shift.01.bam $(OUTDIR)/shift.ref.fa $(OUTDIR)/shift.ref.fa.fai $(OUTDIR)/delly
	$(OUTDIR)/delly call -g $(OUTDIR)/shift.ref.fa -o $@ $(OUTDIR)/shift.01.bam
	

$(OUTDIR)/delly :
	mkdir -p $(dir $@)
	wget -O "$@" "https://github.com/dellytools/delly/releases/download/v1.3.1/delly_v1.3.1_linux_x86_64bit"
	chmod +x $@

$(OUTDIR)/shift.01.bam: $(BAM) $(OUTDIR)/locus.bed
	samtools view -M -L  $(OUTDIR)/locus.bed -O BAM $(BAM)  |\
	java -jar $${JVARKIT_DIST}/jvarkit.jar shiftbam  --bed  $(OUTDIR)/locus.bed |\
	samtools sort -T $(OUTDIR)/tmp --write-index -O BAM -o $@

$(OUTDIR)/shift.ref.dict: $(OUTDIR)/shift.ref.fa
	samtools dict $< > $@

$(OUTDIR)/shift.ref.fa.fai: $(OUTDIR)/shift.ref.fa
	samtools faidx $<

$(OUTDIR)/shift.ref.fa : $(REF)
	mkdir -p $(dir $@)
	samtools faidx $< $(INTERVAL) > $@

$(OUTDIR)/locus.bed:
	mkdir -p $(dir $@)
	echo '$(INTERVAL)' | awk -F '[:-]' '{printf("%s\t%d\t%s\n",$$1,int($$2)-1,$$3);}' > $@

```



END_DOC
*/
@Program(name="shiftbam",
description="shit all coordinates of a bam",
keywords={"bam","sam","bed"},
creationDate="20241001",
modificationDate="20241001",
jvarkit_amalgamion = true
)
public class ShiftBam extends Launcher {
	private static final Logger LOG=Logger.of(ShiftBam.class); 
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"--bed"},description="Bed Path. Extract Reads overlaping this BED. Or use -R2.")
	private Path bedPath = null;
	@Parameter(names={"-R","--reference","--source-reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path sourceRefPath = null;
	@Parameter(names={"-R2","--destination-reference"},description="Original fasta reference. We shift the bam back to this reference. Required without --bed")
	private Path destinationRefPath = null;
	@Parameter(names={"--stringency"},description="Validation Stringency")
	private ValidationStringency validationStringency = ValidationStringency.LENIENT;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	
	
	private static class TidPos {
		String contig = SAMRecord.NO_ALIGNMENT_REFERENCE_NAME;
		int pos = SAMRecord.NO_ALIGNMENT_START;
	}
	
	private static Locatable name2interval(final String s) {
		final int colon=s.lastIndexOf(':');
		if(colon==-1) throw new IllegalArgumentException("cannot find ':' in "+s);
		int hyphen=s.lastIndexOf('-');
		if(hyphen==-1 || hyphen < colon) throw new IllegalArgumentException("cannot find '-' in "+s);
		final String contig = s.substring(0,colon);
		final int start = Integer.parseInt(s.substring(colon+1,hyphen));
		final int end = Integer.parseInt(s.substring(hyphen+1));
		if(start>end) throw new IllegalArgumentException("start> end in "+s);
		return new SimpleInterval(contig, start, end);
		}
	
	
	private static TidPos convert(final Map<Integer,Locatable> tid2interval,int tid,int pos) {
		final TidPos p = new TidPos();
		if(tid==SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) return p;
		final Locatable loc = tid2interval.get(tid);
		if(loc==null) throw new IllegalArgumentException();
		p.contig = loc.getContig();
		p.pos = pos + loc.getStart() -1;
		return p;
		}
	@Override
	public int doWork(final List<String> args) {
		if(this.bedPath==null) {
			return doWorkWithoutBed(args);
			}
		else
			{
			return doWorkWithBed(args);
			}
		}
	
	
	private int doWorkWithoutBed(final List<String> args) {
		try {
			if(this.destinationRefPath==null) {
				LOG.error("destination reference path is required.");
				return -1;
				}
			
			final SAMSequenceDictionary destDict = new SequenceDictionaryExtractor().
					extractRequiredDictionary(this.destinationRefPath);
			
			
			final String input = this.oneFileOrNull(args);
			
			final SamReaderFactory srf = SamReaderFactory.
					make().
					validationStringency(this.validationStringency).
					referenceSequence(this.sourceRefPath);
			final SamInputResource samInputResource = input==null?SamInputResource.of(stdin()):SamInputResource.of(input);
			try(SamReader sr=srf.open(samInputResource)) {
				final SAMFileHeader headerIn = sr.getFileHeader();
				final SAMSequenceDictionary dict0 =  SequenceDictionaryUtils.extractRequired(headerIn);
				if(SequenceUtil.areSequenceDictionariesEqual(dict0, destDict)) {
					LOG.error("dictionaries for input BAM and R2 are the same");
					return -1;
					}
				
				if(sourceRefPath!=null) {
					final SAMSequenceDictionary sourceDict = new SequenceDictionaryExtractor().
							extractRequiredDictionary(this.sourceRefPath);
					
					throw new JvarkitException.DictionariesAreNotTheSame(dict0, sourceDict);
					}
				final Map<Integer,Locatable> tid2interval = new HashMap<>(dict0.size());
				final IntervalTreeMap<Locatable> intervalTreeMap = new IntervalTreeMap<>();
				for(SAMSequenceRecord ssr: dict0.getSequences()) {
					final Locatable loc;
					try {
						loc = name2interval(ssr.getContig());
						}
					catch(Throwable err) {
						throw new IllegalArgumentException("In input SAM :"+ssr.getContig()+" is not a valid contig. expected chrom:start-end", err);
						}
					if(destDict.getSequence(loc.getContig())==null) {
						throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), destDict);
						}
					
					if(intervalTreeMap.containsOverlapping(loc)) {
						throw new IllegalArgumentException("intervals in dictionary overlap and it's forbiden. Intervals should have been merged. see"+loc);
						}
					intervalTreeMap.put(new Interval(loc), loc);
					tid2interval.put(ssr.getSequenceIndex(), loc);
					}
				
				final SAMFileHeader headerOut = headerIn.clone();
				headerOut.setSortOrder(SAMFileHeader.SortOrder.unsorted);
				headerOut.setSequenceDictionary(destDict);
				JVarkitVersion.getInstance().addMetaData(this, headerOut);
								
				try(CloseableIterator<SAMRecord> iter = sr.iterator()) {
					try(SAMFileWriter w= writingBamArgs.setReferencePath(this.destinationRefPath).openSamWriter(this.outputFile, headerOut, true)) {
						while(iter.hasNext()) {
							final SAMRecord rec0 = iter.next();
							final TidPos p1 = convert(tid2interval, rec0.getReferenceIndex(), rec0.getAlignmentStart());
							final TidPos p2 = convert(tid2interval, rec0.getMateReferenceIndex(), rec0.getMateAlignmentStart());
							
							
							final SAMRecord rec1 = rec0.deepCopy();
							rec1.setHeader(headerOut);
							rec1.setReferenceName(p1.contig);
							rec1.setAlignmentStart(p1.pos);
			            	rec1.setMateReferenceName(p2.contig);
			            	rec1.setMateAlignmentStart(p2.pos);
				            
			            	
			            	
			            	// suppl alignments
			                final Object saValue = rec0.getAttribute(SAMTag.SA);
			                if (saValue != null) {
			                	final String saStr = String.class.cast(saValue);
			                	final String[] semiColonStrs =CharSplitter.SEMICOLON.split(saStr);
			                	final List<String> saList = new ArrayList<>(semiColonStrs.length);
			                	for(int i=0;i< semiColonStrs.length;i++) {
			                		final String semiColonStr=semiColonStrs[i];
			                		if(semiColonStr.isEmpty()) continue;
			                		final String commaStrs[] = CharSplitter.COMMA.split(semiColonStr);
			                		
			                		final int saTid = dict0.getSequenceIndex(commaStrs[0]);
			                        if (saTid == -1)  throw new JvarkitException.ContigNotFoundInDictionary(semiColonStr,dict0);
			                        final int  saStart = Integer.parseInt(commaStrs[1]);
			                        final TidPos p3 = convert(tid2interval,saTid,saStart);
			                        commaStrs[0] = p3.contig;
			                        commaStrs[1] = String.valueOf(p3.pos);
			                        saList.add(String.join(",",commaStrs));
			                		}
			                	rec1.setAttribute(SAMTag.SA, String.join(";", saList));
			                	}
							w.addAlignment(rec1);
							}
						}
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	private static int convert2pos(Locatable region,int pos) {
		if(pos< region.getStart() ) throw new IllegalArgumentException("pos < interval.start interval is "+region.toString()+" pos="+pos);
		if(pos> region.getEnd()) throw new IllegalArgumentException("pos > interval.end interval is "+region.toString()+" pos="+pos);
		return 1 + pos - region.getStart();
		}
	
	private int doWorkWithBed(final List<String> args) {
		try {		
			final String input = this.oneFileOrNull(args);
			
			final SamReaderFactory srf = SamReaderFactory.
					make().
					validationStringency(this.validationStringency).
					referenceSequence(this.sourceRefPath);
			final SamInputResource samInputResource = input==null?SamInputResource.of(stdin()):SamInputResource.of(input);
			try(SamReader sr=srf.open(samInputResource)) {
				final SAMFileHeader headerIn = sr.getFileHeader();
				final SAMSequenceDictionary dict0 =  SequenceDictionaryUtils.extractRequired(headerIn);
				
				
				final IntervalTreeMap<Interval> intervalTreeMap=new IntervalTreeMap<>();
				try(BedLineReader rd = new BedLineReader(this.bedPath)) {
					while(rd.hasNext()) {
						final BedLine bedrec = rd.next();
						final Interval loc= new Interval(bedrec.getContig(),bedrec.getStart(),bedrec.getEnd(),false,bedrec.getContig()+":"+bedrec.getStart()+"-"+bedrec.getEnd());
						if(dict0.getSequence(loc.getContig())==null) {
							throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), dict0);
							}
						if(intervalTreeMap.containsOverlapping(loc)) {
							throw new IllegalArgumentException("intervals in "+this.bedPath+" overlap and it's forbiden. Intervals should have been merged. see"+loc);
							}
						intervalTreeMap.put(loc, loc);
						}
					}

				final SAMSequenceDictionary destDict = new SAMSequenceDictionary(
					intervalTreeMap.values().
						stream().
						sorted(new ContigDictComparator(dict0).createLocatableComparator()).
						map(LOC->new SAMSequenceRecord(LOC.getName(), LOC.getLengthOnReference())).
						collect(Collectors.toList())
					);
				
				if(this.sourceRefPath!=null) {
					final SAMSequenceDictionary sourceDict = new SequenceDictionaryExtractor().
							extractRequiredDictionary(this.sourceRefPath);
					
					SequenceUtil.assertSequenceDictionariesEqual(dict0, sourceDict);
					}
				
				final SAMFileHeader headerOut = headerIn.clone();
				headerOut.setSortOrder(SAMFileHeader.SortOrder.unsorted);
				headerOut.setSequenceDictionary(destDict);
				JVarkitVersion.getInstance().addMetaData(this, headerOut);
								
				try(CloseableIterator<SAMRecord> iter = sr.iterator()) {
					try(SAMFileWriter w= writingBamArgs.
								setReferencePath(this.destinationRefPath).
								openSamWriter(this.outputFile, headerOut, true))
						{
						while(iter.hasNext()) {
							final SAMRecord rec0 = iter.next();
							if(rec0.getReadUnmappedFlag()) continue;
							
							final Interval loc0 = intervalTreeMap.getOverlapping(rec0).
									stream().
									filter(RGN->CoordMath.encloses(RGN.getStart(), RGN.getEnd(), rec0.getStart(), rec0.getEnd())).
									findFirst().
									orElse(null);
							if(loc0==null) {
							    continue;
								}
							
							
								
							final SAMRecord rec1 = rec0.deepCopy();
							
							
							
							rec1.setHeader(headerOut);
			            	rec1.setReferenceIndex(destDict.getSequenceIndex(loc0.getName()));
							rec1.setReferenceName(loc0.getName());
							rec1.setAlignmentStart(convert2pos(loc0,rec0.getAlignmentStart()));
							
							if(rec0.getReadPairedFlag()) {
								if(rec0.getMateReferenceIndex()!=SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX/* no !rec0.getMateUnmappedFlag() */) {
									final int mateStart = rec0.getMateAlignmentStart();
									final int mateEnd= rec0.hasAttribute(SAMTag.MC)?SAMUtils.getMateAlignmentEnd(rec0):mateStart;
									final Locatable mateRgn = new SimpleInterval(rec0.getMateReferenceName(),mateStart,mateEnd);
									final Interval loc1= intervalTreeMap.getOverlapping(mateRgn).
											stream().
											filter(RGN->CoordMath.encloses(RGN.getStart(), RGN.getEnd(), mateRgn.getStart(), mateRgn.getEnd())).
											findFirst().
											orElse(null);
									if(loc1==null) {
										rec1.setMateReferenceIndex(SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
										rec1.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
										rec1.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
										}
									else
										{
						            	rec1.setMateReferenceName(loc1.getName());
						            	rec1.setMateReferenceIndex(destDict.getSequenceIndex(loc1.getName()));
						            	rec1.setMateAlignmentStart(convert2pos(loc1,rec0.getMateAlignmentStart()));
										}
									if(loc1==null || !loc1.equals(loc0)) {
										rec1.setProperPairFlag(false);
										rec1.setInferredInsertSize(0);
										}
									
									}
								}
			            	// suppl alignments
			                final Object saValue = rec0.getAttribute(SAMTag.SA);
			                if (saValue != null) {
			                	final String saStr = String.class.cast(saValue);
			                	final String[] semiColonStrs =CharSplitter.SEMICOLON.split(saStr);
			                	final List<String> saList = new ArrayList<>(semiColonStrs.length);
			                	for(int i=0;i< semiColonStrs.length;i++) {
			                		final String semiColonStr=semiColonStrs[i];
			                		if(semiColonStr.isEmpty()) continue;
			                		final String commaStrs[] = CharSplitter.COMMA.split(semiColonStr);
			                		final int saTid = dict0.getSequenceIndex(commaStrs[0]);
			                        if (saTid == -1)  throw new JvarkitException.ContigNotFoundInDictionary(semiColonStr,dict0);
			                        final int  saStart = Integer.parseInt(commaStrs[1]);
			                        final Locatable saLoc = new SimplePosition(commaStrs[0],saStart);
			                        final Interval loc3= intervalTreeMap.getOverlapping(saLoc).
			                        		stream().
			                        		filter(RGN->CoordMath.encloses(RGN.getStart(), RGN.getEnd(), saLoc.getStart(), saLoc.getEnd())).
			                        		findFirst().
			                        		orElse(null);
			                        if(loc3==null) continue;
			                        commaStrs[0] = loc3.getName();
			                        commaStrs[1] = String.valueOf(convert2pos(loc3, saStart));
			                        saList.add(String.join(",",commaStrs));
			                		}
			                	if(saList.isEmpty()) {
			                		rec1.setAttribute(SAMTag.SA, null);
			                		}
			                	else
			                		{
			                		rec1.setAttribute(SAMTag.SA, String.join(";", saList));
			                		}
			                	}
							w.addAlignment(rec1);
							}
						}
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(final String[] args) {
		new ShiftBam().instanceMainWithExit(args);
	}

}
