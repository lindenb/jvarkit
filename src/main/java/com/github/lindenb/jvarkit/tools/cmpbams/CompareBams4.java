/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.cmpbams;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.Counter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ReadNameSortMethod;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.Log.LogLevel;
import htsjdk.samtools.util.PeekableIterator;
/*
BEGIN_DOC

## Example
The following Makefile compare the bam for hg19 and hg38 on chr22 and 21

```Makefile

include ../../config/config.mk

CHROMS=21 22
OUTDIR=tmp

define run
${OUTDIR}/$(1).bam : ${OUTDIR}/$(1).fa.bwt R1.fastq.gz R1.fastq.gz
	${bwa.exe} mem -M   -R '@RG\tID:SAMPLE\tLB:SAMPLE\tSM:SAMPLE\tPL:illumina\tCN:Nantes' ${OUTDIR}/$(1).fa $$(word 2,$$^) $$(word 3,$$^) |  ${samtools.exe} view -b -u -S -F4 - | ${samtools.exe} sort -n -o $$@ -T ${OUTDIR}/$(1)_tmp -

${OUTDIR}/$(1).fa.bwt : ${OUTDIR}/$(1).fa
	${bwa.exe} index $$<

${OUTDIR}/$(1).dict : ${OUTDIR}/$(1).fa
	${java.exe} -jar $(picard.jar) CreateSequenceDictionary R=$$< O=$$@

${OUTDIR}/$(1).fa.fai : ${OUTDIR}/$(1).fa
	${samtools.exe} faidx $$<
	
${OUTDIR}/$(1).fa : 
	mkdir -p $$(dir $$@) && rm -f $$@
	$$(foreach C,${CHROMS}, curl "http://hgdownload.cse.ucsc.edu/goldenPath/$(1)/chromosomes/chr$${C}.fa.gz" | gunzip -c >> $$@;)

endef

all:  ${OUTDIR}/diff.txt 

${OUTDIR}/diff.txt : ${OUTDIR}/hg19.bam ${OUTDIR}/hg38.bam ${OUTDIR}/hg19ToHg38.over.chain ${OUTDIR}/hg19.dict ${OUTDIR}/hg38.dict
	mkdir -p $(dir $@) &&  $(call run_jvarkit,cmpbams4) --novalidchain -st -c $(word 3,$^) $(word 1,$^)  $(word 2,$^) > $@

${OUTDIR}/hg19ToHg38.over.chain :
	mkdir -p $(dir $@) && curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" | gunzip -c > $@

$(eval $(call run,hg19))
$(eval $(call run,hg38))

```

#### Output

```
onlyIn	liftover	compareContig	shift	diffCigarOperations	diffNM	diffFlags	diffChroms	Count
BOTH	SameChrom	DiscordantContig	.	-1	0	147/163	chr22/chr21	2
BOTH	SameChrom	SameContig	Gt100	3	15	83/83	chr22/chr22	1
BOTH	SameChrom	DiscordantContig	.	3	5	147/129	chr21/chr22	1
BOTH	SameChrom	SameContig	Gt100	-1	1	163/163	chr22/chr22	22
BOTH	SameChrom	DiscordantContig	.	0	1	83/99	chr21/chr22	32
BOTH	SameChrom	SameContig	Gt100	0	2	99/99	chr21/chr21	22
BOTH	SameChrom	DiscordantContig	.	2	6	81/65	chr22/chr21	1
BOTH	SameChrom	SameContig	Gt100	0	0	185/137	chr22/chr22	20
BOTH	SameChrom	SameContig	Zero	0	0	177/177	chr22/chr22	1417
(...)
```


END_DOC
 */
@Program(name="cmpbams4",
	description="Compare two BAM files. Print a tab-delimited report",
	keywords={"sam","bam","compare"})
public class CompareBams4  extends Launcher
	{
	private static final Logger LOG = Logger.build(CompareBams4.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-c","--chain"},description="Lift Over file from bam1 to bam2. Optional")
	private File chainFile = null;

	@Parameter(names={"-m","--mismatch"},description="Default Lift Over mismatch. negative=use default")
	private double liftOverMismatch = LiftOver.DEFAULT_LIFTOVER_MINMATCH ;

	@Parameter(names={"-novalidchain","--novalidchain"},description="Disable Lift Over chain validation")
	private boolean disableChainValidation = false;

	@Parameter(names={"-sortmethod","--sortmethod"},description="[20171110]"+ReadNameSortMethod.DESCRIPTION)
	private ReadNameSortMethod sortMethod = ReadNameSortMethod.picard;

	private LiftOver liftOver=null;
	private enum OnlyIn { ONLY_IN_FIRST,ONLY_IN_SECOND,BOTH};
	private enum LiftOverStatus {NoLiftOver,SourceUnmapped,DestNotInDict,LiftOverFailed,DiscordantChrom,SameChrom};
	private enum CompareContig {BothUnmapped,GainMapping,LostMapping,SameContig,DiscordantContig};
	private enum Shift {Zero,Le10,Le20,Le100,Gt100};
	private enum Flag { Discordant,Same};
		
	private static int hash(final int prev,final Object o) {
		return prev*31 +(o==null?0:o.hashCode());
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	private static int cmp(final Comparable a,final Comparable b) {
		if(a==null && b!=null) return -1;
		if(a!=null && b==null) return 1;
		if(a==null && b==null) return 0;
		return a.compareTo(b);
	}
	
	private static void str(final StringBuilder sb,final Object o) {
		sb.append(o==null?".":o.toString()).append("\t");
	}
	
	private static class Couple<T extends Comparable<T>> implements Comparable<Couple<T>>
	{
		final T v1;
		final T v2;
		public Couple(final T v1,final T v2) {
			this.v1 = v1;
			this.v2 = v2;
		}
		
		@Override
		public int compareTo(final Couple<T> o) {
			int i=v1.compareTo(o.v1);if(i!=0) return i;
			i=v2.compareTo(o.v2);if(i!=0) return i;
			return 0;
			}
		
		public int hashCode() { return v1.hashCode()*31+v2.hashCode();}
		@SuppressWarnings("unchecked")
		@Override
		public boolean equals(final Object obj) {
			if( obj==this) return true;
			if(obj == null || !(obj instanceof Couple)) return false;
			return this.compareTo(Couple.class.cast(obj)) == 0;
		}

		@Override
		public String toString() {
			return String.valueOf(v1)+"/"+v2;
			}
		
		
	}

	
	private static final class Diff
		implements Comparable<Diff>
		{
		OnlyIn onlyIn = null;
		LiftOverStatus liftover = null;
		CompareContig compareContig = null;
		Shift shift = null;
		Integer diffCigarOperations = null;
		Integer diffNM = null;
		Couple<Integer> diffFlags = null;
		Couple<String> diffChrom = null;
		Flag diffFlag;
		
		@Override
		public int hashCode() {
			int h = 0;
			h = hash(h,onlyIn);
			h = hash(h,liftover);
			h = hash(h,compareContig);
			h = hash(h,shift);
			h = hash(h,diffCigarOperations);
			h = hash(h,diffNM);
			h = hash(h,diffFlags);
			h = hash(h,diffChrom);
			h = hash(h,diffFlag);
			return h;
		}
		
		@Override
		public int compareTo(final Diff o) {
			if( o==this) return 0;
			int i;
			i = cmp(this.onlyIn,o.onlyIn); if(i!=0) return i;
			i = cmp(this.liftover,o.liftover); if(i!=0) return i;
			i = cmp(this.compareContig,o.compareContig); if(i!=0) return i;
			i = cmp(this.shift,o.shift); if(i!=0) return i;
			i = cmp(this.diffCigarOperations,o.diffCigarOperations); if(i!=0) return i;
			i = cmp(this.diffNM,o.diffNM); if(i!=0) return i;
			i = cmp(this.diffFlags,o.diffFlags); if(i!=0) return i;
			i = cmp(this.diffChrom,o.diffChrom); if(i!=0) return i;
			i = cmp(this.diffFlag,o.diffFlag); if(i!=0) return i;
			return 0;
			}
		
		@Override
		public boolean equals(final Object obj) {
			if( obj==this) return true;
			if(obj == null || !(obj instanceof Diff)) return false;
			return this.compareTo(Diff.class.cast(obj)) == 0;
			}
		
		@Override
		public String toString() {
			final StringBuilder sb = new StringBuilder();
			str(sb,this.onlyIn);
			str(sb,this.liftover);
			str(sb,this.compareContig);
			str(sb,this.shift);
			str(sb,this.diffCigarOperations);
			str(sb,this.diffNM);
			str(sb,this.diffFlags);
			str(sb,this.diffChrom);
			str(sb,this.diffFlag);
			return sb.toString();
			}
		}
	
	private static Interval interval(final SAMRecord rec) {
		if(rec.getReadUnmappedFlag()) return null;
		return new Interval(
			rec.getReferenceName(),
			rec.getAlignmentStart(),
			rec.getAlignmentEnd(),
			rec.getReadNegativeStrandFlag(),
			null
			);
	}
	
	/** used by Comm Bams */
	static int side(final SAMRecord rec) {
		if(rec.getReadPairedFlag()) {
			if(rec.getFirstOfPairFlag()) return 1;
			if(rec.getSecondOfPairFlag()) return 2;
			throw new RuntimeException("Side for flag "+rec.getReadName()+":"+rec.getFlags()+"?");
		}
		else {
			return 0;
		}
	}
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter out=null;
		final SamReader samFileReaders[]=new SamReader[]{null,null};
		final SAMSequenceDictionary dicts[]=new SAMSequenceDictionary[]{null,null};
		@SuppressWarnings("unchecked")
		final PeekableIterator<SAMRecord> iters[]=new PeekableIterator[]{null,null};
		final List<List<SAMRecord>> recordLists = Arrays.asList(new ArrayList<SAMRecord>(),new ArrayList<SAMRecord>());
		final Counter<Diff> diffs = new Counter<>();
		final short NM_TAG = SAMTag.NM.getBinaryTag();
		if(args.size() !=2)
			{
			LOG.error("Expected two and only two bams please, but got "+args.size());
			return -1;
			}
		
		try
			{
			final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			
			for(int i=0;i< args.size() && i< samFileReaders.length;++i)
				{
				final File samFile = new File(args.get(i));
				LOG.info("opening "+samFile);
				samFileReaders[i]=srf.open(samFile);
				final SAMFileHeader header = samFileReaders[i].getFileHeader();
				if(header.getSortOrder()!=SAMFileHeader.SortOrder.queryname) {
					LOG.error("Expected "+samFile+" to be sorted on "+SAMFileHeader.SortOrder.queryname+" but got "+header.getSortOrder());
					return -1;	
					}
				dicts[i] = header.getSequenceDictionary();
				if( dicts[i]==null || dicts[i].isEmpty()) {
					LOG.error("In "+samFile+": No SAMSequenceDictionary in header.");
					return -1;
					}
				iters[i] = new PeekableIterator<>(samFileReaders[i].iterator());
				}
			
			if(this.chainFile!=null) {
				LOG.info("opening chain file "+this.chainFile+".");
				this.liftOver =new LiftOver(this.chainFile);
				this.liftOver.setLiftOverMinMatch(this.liftOverMismatch<=0?LiftOver.DEFAULT_LIFTOVER_MINMATCH:this.liftOverMismatch);

				if(!this.disableChainValidation) {
				this.liftOver.validateToSequences(dicts[1]);
				}
			}
			
			final Comparator<SAMRecord> comparator = this.sortMethod.get();
			
			for(;;) {
				for(int i=0;i< 2;++i) {
					if(recordLists.get(i).isEmpty())
						{
						while(iters[i].hasNext()) {
							final SAMRecord rec = iters[i].peek();
							if(rec.isSecondaryOrSupplementary()) {
								iters[i].next();
								continue;
								}
							if(!recordLists.get(i).isEmpty() && comparator.compare(recordLists.get(i).get(0),rec)>0)
								{
								throw new SAMException("Something is wrong in sort order of "+args.get(i)+" : got\n\t"
										+rec+" "+side(rec)+"\nafter\n\t"+ recordLists.get(i).get(0)+" "+side(recordLists.get(i).get(0))+"\nSee also option (samtools querysort)"
										);
								}
							else if( recordLists.get(i).isEmpty() ||
								comparator.compare(recordLists.get(i).get(0),rec)==0)
								{
								recordLists.get(i).add(iters[i].next());
								}
							else
								{	
								break;
								}
							}
						}
				}
				if(recordLists.get(0).isEmpty() && recordLists.get(1).isEmpty() ) break;
				
				if(recordLists.get(0).size()>1) {
					LOG.warn("size>2 for 1/2:"+recordLists.get(0).get(0).getReadName());
					for(final SAMRecord sr:recordLists.get(0))
						{	
						LOG.warn(">> "+sr+" flags:"+sr.getFlags()+" pos:"+sr.getReferenceName()+":"+sr.getStart());
						}
					}
				if(recordLists.get(1).size()>1) {
					LOG.warn("size>2 for 2/2:"+recordLists.get(1).get(0).getReadName());
					for(final SAMRecord sr:recordLists.get(1))
						{	
						LOG.warn(">> "+sr+" flags:"+sr.getFlags()+" pos:"+sr.getReferenceName()+":"+sr.getStart());
						}

				}
				
				final Diff diff = new Diff();
				
								
				if((recordLists.get(0).isEmpty() && !recordLists.get(1).isEmpty()) ||
				   (!recordLists.get(0).isEmpty() && !recordLists.get(1).isEmpty() && comparator.compare(recordLists.get(0).get(0), recordLists.get(1).get(0))>0)
					) {
					diff.onlyIn = OnlyIn.ONLY_IN_SECOND;
					recordLists.get(1).clear();
					}
				else if((!recordLists.get(0).isEmpty() && recordLists.get(1).isEmpty()) ||
						(!recordLists.get(0).isEmpty() && !recordLists.get(1).isEmpty() && comparator.compare(recordLists.get(0).get(0), recordLists.get(1).get(0))<0))
					{
					diff.onlyIn = OnlyIn.ONLY_IN_FIRST;
					recordLists.get(0).clear();
					}
				else
				{
					final SAMRecord rec0 = recordLists.get(0).get(0);
					final SAMRecord rec1 = recordLists.get(1).get(0);
					final Interval rgn0a = interval(rec0);
					final Interval rgn0b;
					
					diff.onlyIn = OnlyIn.BOTH;
					diff.diffFlags = new Couple<Integer>(rec0.getFlags(),rec1.getFlags());
					diff.diffFlag = (rec0.getFlags()==rec1.getFlags()?Flag.Same:Flag.Discordant);

					
					if(this.liftOver==null) {
						rgn0b = rgn0a;
						diff.liftover = LiftOverStatus.NoLiftOver;
						}
					else if(rgn0a==null) {
						diff.liftover = LiftOverStatus.SourceUnmapped;
						rgn0b = rgn0a;
						}
					else {
						rgn0b = this.liftOver.liftOver(rgn0a);
						if(rgn0b==null) {
							diff.liftover = LiftOverStatus.LiftOverFailed;
							}
						else if(dicts[1].getSequence(rgn0b.getContig())==null)
							{
							diff.liftover = LiftOverStatus.DestNotInDict;
							}
						else if(rgn0a.getContig().equals(rgn0b.getContig())) {
							diff.liftover = LiftOverStatus.SameChrom;
							}
						else
							{
							diff.liftover = LiftOverStatus.DiscordantChrom;
							}
						}
					
					final Interval rgn1 = interval(rec1);
					if(rgn0b==null && rgn1==null) {
						diff.compareContig = CompareContig.BothUnmapped;
					} else if(rgn0b==null && rgn1!=null) {
						diff.compareContig = CompareContig.GainMapping;
						diff.diffChrom = new Couple<String>("*",rgn1.getContig());

					} else if(rgn0b!=null && rgn1==null) {
						diff.compareContig = CompareContig.LostMapping;
						diff.diffChrom = new Couple<String>(rgn0b.getContig(),"*");
					} else if(rgn0b.getContig().equals(rgn1.getContig()))
					{
						diff.compareContig = CompareContig.SameContig;
						diff.diffChrom = new Couple<String>(rgn0b.getContig(),rgn1.getContig());
						final int shift = Math.abs(rgn0b.getStart() - rgn1.getStart());
						if(shift == 0 ) {
							diff.shift = Shift.Zero;
						}
						else if(shift <=10) {
							diff.shift = Shift.Le10;
							
						}else if(shift <=20) {
							diff.shift = Shift.Le20;
						}
						else if( shift<= 100)
						{
							diff.shift = Shift.Le100;
						} else {
							diff.shift = Shift.Gt100;
						}
						
					} else
					{
						
						if(dicts[1].getSequence(rgn0b.getContig())==null)
							{
							diff.diffChrom = new Couple<String>("?",rgn1.getContig());
							}
						else
							{
							diff.diffChrom = new Couple<String>(rgn0b.getContig(),rgn1.getContig());
							}
						diff.compareContig = CompareContig.DiscordantContig;
					}
					
					if(rgn0b!=null && rgn1!=null) {
						diff.diffCigarOperations = (
								rec0.getCigar().numCigarElements() - rec1.getCigar().numCigarElements()
								);
						final Object nm0 = rec0.getAttribute(NM_TAG);
						final Object nm1 = rec1.getAttribute(NM_TAG);
						if(nm0!=null && nm1!=null) {
							diff.diffNM = (Number.class.cast(nm0).intValue()-Number.class.cast(nm1).intValue());
						}
					}
					recordLists.get(1).clear();
					recordLists.get(0).clear();
					}
				
				diffs.incr(diff);
			}
			
			
			LOG.info("done");
			
			final StringBuilder sb = new StringBuilder();
			str(sb,"onlyIn");
			str(sb,"liftover");
			str(sb,"compareContig");
			str(sb,"shift");
			str(sb,"diffCigarOperations");
			str(sb,"diffNM");
			str(sb,"diffFlags");
			str(sb,"diffChroms");
			str(sb,"diffFlag");
			sb.append("Count");
			out = super.openFileOrStdoutAsPrintWriter(outputFile);
			out.println(sb);
			for(final Diff key:diffs.keySet()) {
				out.print(key.toString());
				out.println(diffs.count(key));
			}
			out.flush();
			out.close();
			
			return RETURN_OK;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			this.liftOver = null;
			for(int i=0;i< samFileReaders.length;++i){
				CloserUtil.close(iters[i]);
				CloserUtil.close(samFileReaders[i]);
				}
			CloserUtil.close(out);
			}
		}
		
	public static void main(String[] args) throws Exception
		{
		 Log.setGlobalLogLevel(LogLevel.ERROR);
		new CompareBams4().instanceMainWithExit(args);
		}
}
