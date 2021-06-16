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
package com.github.lindenb.jvarkit.tools.pcr;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassBamLauncher;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.IntervalListProvider;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;


/**
BEGIN_DOC


## See also

* (2021) `samtools ampliconclip` – clip reads using a BED file  http://www.htslib.org/doc/samtools-ampliconclip.html


## Input/Output

input is a bam

output is a bam


name of the BED record is appended to the original  read name/

unmapped reads, reads without cigar or reads that don't overlap any BED record are discarded

MAPQ is set to 255

reads are converted to singled end

optional args are not filled

bounding bases with cigar string without cigar operator M/X/= are discarded.



## Example

```
$ cat jeter.bed
RF01	10	15
RF01	20	25
RF01	30	35
```

```
$ java -jar dist/bamslicebed.jar -B jeter.bed ./src/test/resources/S1.bam |samtools sort -T tmp -o jeter.bam -
$ samtools view jeter.bam 
RF01_1_483_2:0:0_3:0:0_41#RF01:11:15	0	RF01	11	255	5M	*	0	0	GCTAT	22222
RF01_8_542_1:0:0_2:0:0_95#RF01:11:15	0	RF01	11	255	5M	*	0	0	GCTAT	22222
RF01_11_507_0:0:0_1:0:0_9e#RF01:11:15	0	RF01	11	255	5M	*	0	0	GCTAT	22222
RF01_12_501_0:0:0_2:0:0_62#RF01:11:15	0	RF01	12	255	4M	*	0	0	CTAT	2222
RF01_1_483_2:0:0_3:0:0_41#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGC	22222
RF01_8_542_1:0:0_2:0:0_95#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGA	22222
RF01_11_507_0:0:0_1:0:0_9e#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGA	22222
RF01_12_501_0:0:0_2:0:0_62#RF01:21:25	0	RF01	21	255	5M	*	0	0	GGGGA	22222
RF01_1_483_2:0:0_3:0:0_41#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_8_542_1:0:0_2:0:0_95#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_11_507_0:0:0_1:0:0_9e#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_12_501_0:0:0_2:0:0_62#RF01:31:35	0	RF01	31	255	5M	*	0	0	AATCT	22222
RF01_27_590_3:0:0_1:0:0_68#RF01:31:35	0	RF01	31	255	5M	*	0	0	CATCT	22222

samtools tview jeter.bam src/test/resources/rotavirus_rf.fa


1         11        21        31    
ggctattaaagctatacaATGGGGAAGTATAATCTA
          .....     .....     .....
          .....     ....C     .....
          .....     .....
          .....     .....
           ....     .....
                              .....
                              .....
                              .....
                              C....
```

END_DOC
*/

@Program(name="bamslicebed",
	description="For @wouter_decoster : slice (long reads) overlapping the records of a BED file",
	keywords={"sam","bam","bed"},
	creationDate="20191030",
	modificationDate="20210615"
	)
public class BamSliceBed extends OnePassBamLauncher {
	private static final Logger LOG = Logger.build(BamSliceBed.class).make();

	@Parameter(names={"-B","--bed","--pcr"},description="Regions containing non-overlapping PCR fragments. "+IntervalListProvider.OPT_DESC,converter=IntervalListProvider.StringConverter.class,required=true)
	private IntervalListProvider intervalListProvider = IntervalListProvider.unspecified();
	@Parameter(names={"--attributes"},description="Leep the following attributes (separated by spaces/comma/semicolon)")
	private String keepAttributesStr = "";
	@Parameter(names={"--clip"},description="Do not remove the bases but soft clip them.")
	private boolean use_clip = false;

	
	
	private static final byte NO_BASE = '\0';
	private static final char NO_QUAL = '\0';
	private SAMProgramRecord spr = null;
	private final IntervalTreeMap<Interval> bedIntervals = new IntervalTreeMap<>();
	private final Set<String> keepAttributes = new HashSet<>();
	
	private static class Base
		{
		byte readbase=NO_BASE;
		char readqual=NO_QUAL;
		int  readpos=-1;
		int  refpos=-1;
		CigarOperator cigaroperator = null;
		}

	@Override
	protected Logger getLogger() {
		return LOG;
		}	
	
	@Override
	protected int beforeSam()
		{
		this.keepAttributes.addAll(
				Arrays.stream(this.keepAttributesStr.split("[,; \t]")).
				filter(S->!StringUtils.isBlank(S)).
				filter(S->!S.equals("RG")).//always done
				collect(Collectors.toSet()));
		
		for(final String att: this.keepAttributes) {
			if(att.length()!=2) {
				LOG.error("attribute "+att+" doesn't have length=2");
				return -1;
			}
		}
		
		return super.beforeSam();
		}
	
	@Override
	protected SAMFileHeader createOutputHeader(SAMFileHeader headerIn) {
		SAMFileHeader header2 =  super.createOutputHeader(headerIn);
		
		final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(headerIn);
		

		this.intervalListProvider.
			dictionary(dict).
			skipUnknownContigs().
			stream().
			map(R->new Interval(R)).
			forEach(R->bedIntervals.put(R,R));
			
		header2.setSortOrder(SortOrder.unsorted);
		
		this.spr = super.createProgramRecord(header2);
		return header2;
		}
	
	@Override
	protected void scanIterator(final SAMFileHeader headerIn,final CloseableIterator<SAMRecord> iter,final SAMFileWriter sfw) {
			final SAMRecordFactory samRecordFactory = new DefaultSAMRecordFactory();
			while(iter.hasNext())
				{
				final SAMRecord rec=  iter.next();
				if(rec.getReadUnmappedFlag()) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				final Collection<Interval> beds = bedIntervals.getOverlapping(rec);
				if(beds.isEmpty()) continue;
				
				final List<Base> align = new ArrayList<>(rec.getLengthOnReference());
				int refpos=rec.getUnclippedStart();
				int readpos=0;
				final byte bases[] = rec.getReadBases();
				if(bases == SAMRecord.NULL_SEQUENCE) continue;
				
				final String quals;
				final boolean qual_is_available;
				if(rec.getBaseQualities() == SAMRecord.NULL_QUALS) {
					quals = StringUtils.repeatCharNTimes('#',bases.length);
					qual_is_available = false;
					}
				else
					{
					quals = rec.getBaseQualityString();
					qual_is_available = true;
					}
				
				
				for(final CigarElement ce:cigar)
					{
					final CigarOperator op = ce.getOperator();
					switch(op)
						{
						case P: break;
						case D: case N:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.cigaroperator = op;
								b.refpos  = refpos;
								align.add(b);
								refpos++;
								}
							break;
							}
						case H:/* ignore hard clip */
							{
							refpos+=ce.getLength();
							break;
							}
						case S:case X:case EQ:case M:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.refpos=refpos;
								b.cigaroperator = op;
								b.readbase = bases[readpos];
								b.readqual = quals.charAt(readpos);
								b.readpos=readpos;
								align.add(b);
								readpos++;
								refpos++;
								}
							break;
							}
						case  I:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.cigaroperator = op;
								b.readbase = bases[readpos];
								b.readqual = quals.charAt(readpos);
								b.readpos =readpos;
								align.add(b);
								readpos++;
								}
							break;
							}
						default: throw new IllegalStateException();
						}
					}
				
				
				for(final Interval bed: beds) {
					final List<Base> copy = new ArrayList<>(align);
					
					/* trim 5 prime */
					int trim=-1;
					for(int x=0;x< copy.size();x++) {
						final Base b = copy.get(x);
						if(b.cigaroperator.isAlignment() && b.refpos>= bed.getStart()) {
							trim=x;
							break;
							}
						}
					if(trim>0) {
						final List<Base> subList = copy.subList(0, trim);
						if(use_clip){
							for(Base b: subList) {
								if(b.readpos!=-1) {
									b.cigaroperator= CigarOperator.SOFT_CLIP;
									}
								}
							subList.removeIf(B->B.readpos==-1);// N D
							}
						else
							{
							subList.clear();
							}
						}
					
					/* trim 3 prime */
					trim=-1;
					for(int x=copy.size()-1;x>=0;x--) {
						final Base b = copy.get(x);
						if(b.cigaroperator.isAlignment() && b.refpos<= bed.getEnd()) {
							trim=x;
							break;
							}
						}
					if(trim!=-1) {
						final List<Base> subList = copy.subList(trim+1,copy.size());
						if(use_clip) {
							for(Base b: subList) {
							if(b.readpos!=-1) {
								b.cigaroperator= CigarOperator.SOFT_CLIP;
								}
							}
							subList.removeIf(B->B.readpos==-1);// N D
							}
						else {
							subList.clear();
							}
						}
					
					if(copy.stream().noneMatch(P->P.cigaroperator.isAlignment())) continue;
					if(copy.isEmpty()) continue;
					int nrefpos = copy.stream().
							filter(B->B.refpos!=-1).
							filter(B->!B.cigaroperator.isClipping()).
							mapToInt(B->B.refpos).
							findFirst().
							orElse(-1);
					final SAMRecord newrec = samRecordFactory.createSAMRecord(sfw.getFileHeader());
					newrec.setReadName(rec.getReadName()+"#"+bed.getContig()+":"+bed.getStart()+":"+bed.getEnd());
					newrec.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY);
					newrec.setAlignmentStart(nrefpos);
					newrec.setReadString(
							copy.stream().
							filter(B->B.readbase!=NO_BASE).
							map(B->String.valueOf((char)B.readbase)).
							collect(Collectors.joining()));
					
					if(qual_is_available) {
						newrec.setBaseQualityString(
								copy.stream().
								filter(B->B.readqual!=NO_QUAL).
								map(B->String.valueOf(B.readqual)).
								collect(Collectors.joining()));
						}
					else
						{
						newrec.setBaseQualities(SAMRecord.NULL_QUALS);
						}
					newrec.setReadUnmappedFlag(false);
					newrec.setReadNegativeStrandFlag(rec.getReadNegativeStrandFlag());
					newrec.setReferenceName(rec.getReferenceName());
					newrec.setCigar(Cigar.fromCigarOperators(copy.stream().map(B->B.cigaroperator).collect(Collectors.toList())));
					
					if(rec.hasAttribute("RG")) {
						newrec.setAttribute("RG", rec.getAttribute("RG"));
					}
					for(final String att:this.keepAttributes) {
						if(rec.hasAttribute(att)) {
							newrec.setAttribute(att, rec.getAttribute(att));
							}
						}
					
					newrec.setAttribute("PG",this.spr.getId());
					
					
					
					sfw.addAlignment(newrec);
					} //end for interval
				}//end while
			}
	
	public static void main(final String[] args) {
		new BamSliceBed().instanceMainWithExit(args);
		}

}
