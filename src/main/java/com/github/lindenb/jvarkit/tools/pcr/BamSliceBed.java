/*
The MIT License (MIT)

Copyright (c) 2018 Pierre Lindenbaum

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

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IntervalTreeMap;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;


/**
BEGIN_DOC

unmapped reads or reads that don't overlap any BED record are discarded



	* see 


END_DOC
*/

@Program(name="bamslicebed",
	description="For @wouter_decoster : slice (long reads) overlapping the records of a BED file",
	keywords={"sam","bam","bed"}
	)
public class BamSliceBed extends Launcher
	{
	private static final Logger LOG = Logger.build(BamSliceBed.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-B","--bed"},description="Bed file used to slice the bam",required=true)
	private File bedFile = null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	
	private static final byte NO_BASE = '\0';
	private static final byte NO_QUAL = '\0';
	
	private static class Base
		{
		byte readbase=NO_BASE;
		byte readqual=NO_QUAL;
		int  readpos=-1;
		int  refpos=-1;
		CigarOperator cigaroperator = null;
		}

	
	
	@Override
	public int doWork(final List<String> args) {
		if(this.bedFile==null)
			{
			LOG.error("undefined bed file ");
			return -1;
			}
		BufferedReader r=null;
		SamReader samReader=null;
		SAMFileWriter sw=null;
		SAMRecordIterator iter = null;
		try {
			final String inputName = oneFileOrNull(args);
			samReader = openSamReader(inputName);
			final SAMFileHeader header = samReader.getFileHeader();
			if(header==null)
				{
				LOG.error("No SAM header in input");
				return -1;
				}
			final SAMSequenceDictionary dict = header.getSequenceDictionary();
			if(dict==null)
				{
				LOG.error(JvarkitException.BamDictionaryMissing.getMessage(String.valueOf(inputName)));
				return -1;
				}
			final IntervalTreeMap<BedLine> bedIntervals = new IntervalTreeMap<>();

			LOG.info("reading bed File "+this.bedFile);
			r= IOUtils.openFileForBufferedReading(this.bedFile);
			String line;
			final BedLineCodec codec = new BedLineCodec();
			while((line=r.readLine())!=null)
				{
				final BedLine bedLine = codec.decode(line);
				if(bedLine==null)
					{
					LOG.warn("Ignoring bed line "+line);
					continue;
					}
				
				if(dict.getSequenceIndex(bedLine.getContig())==-1) {
					throw new JvarkitException.ContigNotFoundInDictionary(bedLine.getContig(), dict);
					}
					
				bedIntervals.put(bedLine.toInterval(),bedLine);
				}
			CloserUtil.close(r);r=null;
			
			
			final SAMFileHeader header2 = header.clone();
			
			header2.addComment(getProgramName()+" "+getVersion()+": Processed with "+getProgramCommandLine());
			header2.setSortOrder(SortOrder.unsorted);
			
			sw = this.writingBamArgs.openSAMFileWriter(outputFile,header2, false);
			final SAMRecordFactory samRecordFactory = new DefaultSAMRecordFactory();
			final ProgressFactory.Watcher<SAMRecord> progress = ProgressFactory.
							newInstance().
							dictionary(dict).
							logger(LOG).
							build();
			iter =  samReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec= progress.apply(iter.next());
				if(rec.getReadUnmappedFlag()) continue;
				final Cigar cigar = rec.getCigar();
				if(cigar==null || cigar.isEmpty()) continue;
				final Collection<BedLine> beds = bedIntervals.getOverlapping(rec);
				if(beds.isEmpty()) continue;
				
				final List<Base> align = new ArrayList<>();
				int refpos=rec.getUnclippedStart();
				int readpos=0;
				final byte bases[] = rec.getReadBases();
				if(bases == SAMRecord.NULL_SEQUENCE) continue;
				final byte quals[] = rec.getBaseQualities();
				if(quals == SAMRecord.NULL_QUALS) continue;
				
				
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
						case S:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.refpos=refpos;
								b.readpos= readpos;
								b.cigaroperator = op;
								align.add(b);
								refpos++;
								readpos++;
								}
							break;
							}
						case H:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.refpos=refpos;
								b.cigaroperator = op;
								align.add(b);
								refpos++;
								}
							break;
							}
						case X:case EQ:case M:
							{
							for(int i=0;i< ce.getLength();++i)
								{
								final Base b=new Base();
								b.refpos=refpos;
								b.cigaroperator = op;
								b.readbase = bases[readpos];
								b.readqual = quals[readpos];
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
								b.readqual = quals[readpos];
								b.readpos =readpos;
								align.add(b);
								readpos++;
								}
							break;
							}
						default: throw new IllegalStateException();
						}
					}
			

				
				
				for(final BedLine bed: beds) {
					final LinkedList<Base> copy = new LinkedList<>(align);
					Predicate<Base> canRemoveBase = B->{
						if(B.refpos==-1) return true;
						if(B.readpos< bed.getStart()) return true;
						if(B.readpos> bed.getEnd()) return true;
						if(!B.cigaroperator.isAlignment())return true;
						return false;
						};
					while(!copy.isEmpty())
						{
						final Base first = copy.getFirst();
						if(!canRemoveBase.test(first)) break;
						copy.removeFirst();
						}
					while(!copy.isEmpty())
						{
						final Base last = copy.getLast();
						if(!canRemoveBase.test(last)) break;
						copy.removeLast();
						}
					if(copy.stream().noneMatch(P->P.cigaroperator.isAlignment())) continue;
					if(copy.isEmpty()) continue;
					int nrefpos = copy.stream().filter(B->B.refpos!=-1).mapToInt(B->B.readpos).findFirst().orElse(-1);
					final SAMRecord newrec = samRecordFactory.createSAMRecord(header2);
					newrec.setReadName(rec.getReadName()+"#"+bed.getContig()+":"+bed.getStart()+":"+bed.getEnd());
					newrec.setMappingQuality(SAMRecord.UNKNOWN_MAPPING_QUALITY);
					newrec.setAlignmentStart(nrefpos);
					newrec.setReadString(
							copy.stream().
							filter(B->B.readbase!=NO_BASE).
							map(B->String.valueOf((char)B.readbase)).
							collect(Collectors.joining()));
					newrec.setBaseQualityString(
							copy.stream().
							filter(B->B.readqual!=NO_QUAL).
							map(B->String.valueOf((char)B.readqual)).
							collect(Collectors.joining()));
					newrec.setReadUnmappedFlag(false);
					newrec.setReadNegativeStrandFlag(rec.getReadNegativeStrandFlag());
					newrec.setReferenceName(rec.getReferenceName());
					newrec.setCigar(Cigar.fromCigarOperators(copy.stream().map(B->B.cigaroperator).collect(Collectors.toList())));
					
					sw.addAlignment(newrec);
					}
				
			
				}
			sw.close();
			sw=null;
			samReader.close();
			samReader=null;
			progress.close();
			return 0;
			}
		catch (final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(samReader);
			}
		}

	
	public static void main(final String[] args) {
		new BamSliceBed().instanceMainWithExit(args);
		}

}
