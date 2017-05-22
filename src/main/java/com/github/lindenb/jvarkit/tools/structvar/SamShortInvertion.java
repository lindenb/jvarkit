/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SequenceUtil;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.semontology.Term;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlign;
import com.github.lindenb.jvarkit.util.picard.OtherCanonicalAlignFactory;
@Program(name="samshortinvert",
	description="Scan short inversions in SAM",
	terms=Term.ID_0000023
	)
public class SamShortInvertion extends Launcher
	{
	private static final Logger LOG = Logger.build(SamShortInvertion.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-m","--maxsize"},description="max size of inversion")
	private int max_size_inversion = 2000 ;

	@Parameter(names={"-c","--mincov"},description="min coverage")
	private int min_coverage = 10 ;

	@SuppressWarnings("unused")
	private class SamRecordPair
		{
		SAMRecord first;
		SAMRecord second;
		SamRecordPair(SAMRecord first,SAMRecord second)
			{
			if(second==null) {
				this.first=first;
			} else if(first==null) {
				this.first=second;
			  } else if(first.getReferenceName().equals(second.getReferenceName())) {
				  if(first.getAlignmentStart()< second.getAlignmentStart()) {
					  this.first=first;
					  this.second=second;
				  } else
				  	{
					  this.first=second;
					  this.second=first;
				  	}
			  } else
			  	{
				  this.first=first;
				  this.second=second;
			  	}
			}
		SamRecordPair(SAMRecord one)
			{
			
			}
		}
	
	
	private static class SamReaderList
		implements Closeable
		{
		private final List<SamReader> samReaders;
		private final SAMSequenceDictionary dict;
		private final String sample ;
		@SuppressWarnings("unused")
		final File bamFile;
		SamReaderList(final File bamFile,int n) {
			final SamReaderFactory srf = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.LENIENT);
			this.bamFile=bamFile;
			SamReader samReader1 = srf.open(bamFile);
			if(!samReader1.hasIndex()) throw new RuntimeIOException("Bam "+bamFile+" is missing an index");
			SAMFileHeader header= samReader1.getFileHeader();
			this.dict = header.getSequenceDictionary();
			if(this.dict==null)  throw new RuntimeIOException("Bam "+bamFile+" is missing a sequence dictionary");
			String sampleName=null;
			for(final SAMReadGroupRecord srg : header.getReadGroups()) {
				if(srg.getSample()==null) throw new RuntimeIOException("SamReadGroup in"+bamFile+" is missing a SAMPLE (SN)");
				final String indi =srg.getSample();
				if(sampleName!=null && !sampleName.equals(indi)) throw new RuntimeIOException("SamReadGroup in"+bamFile+" : more than one sample defined");
				sampleName=indi;
			}
			if(sampleName==null)  throw new RuntimeIOException("no SamReadGroup in "+bamFile+" with a SAMPLE (SN)");
			this.sample = sampleName;
			List<SamReader> readers = new ArrayList<>();
			readers.add(samReader1);
			while(n>1) {
				readers.add(srf.open(bamFile));
				--n;
				}
			this.samReaders = Collections.unmodifiableList(readers);
			}
		
		public SamReader get(int index)  {
			return this.samReaders.get(index);
		}
		
		@Override
		public void close() throws IOException {
			for(final SamReader r:this.samReaders) r.close();
			}
	}

	
	private void challenge(final SAMRecord rec1, final SAMRecord rec2) {
		if(rec2!=null) {
			if(rec2.getReferenceName().equals(rec1.getReferenceName()) &&
				rec2.getAlignmentStart()<rec1.getAlignmentStart()) {
				challenge(rec2,rec1);
				return;
			}
			if(rec2.getReferenceName().compareTo(rec1.getReferenceName())<0) {
				challenge(rec2,rec1);
				return;
				}
			}
		
		
		}
	
	@Override
	public int doWork(List<String> input) {
		SamReader r=null;
		PrintStream out = null;
		//SAMFileWriter w=null;
		try
			{
			final List<SamReaderList> samReaders = new ArrayList<>();
			final Set<String> args = IOUtils.unrollFiles(input);
			if(args.isEmpty()) {
				LOG.error("No input file");
				return -1;
				}
			SAMSequenceDictionary dict=null;
			for(final String bam: args ) {
				final SamReaderList samReader = new SamReaderList(new File(bam),2);
				for(int i=0;i< samReaders.size();++i )
					{
					if(samReaders.get(i).sample.equals(samReader.sample))
						{
						samReader.close();
						LOG.error("Sample defined in two bams "+samReader.sample);
						return -1;
						}
					if(dict==null) {
						dict = samReader.dict;
					} else if(!SequenceUtil.areSequenceDictionariesEqual(dict, samReader.dict))
						{
						samReader.close();
						LOG.error("bam contains two sequence dict.");
						return -1;
						}
					}
				samReaders.add(samReader);
				}
			
			final Predicate<SAMRecord> samRecordFilter = new Predicate<SAMRecord>() {
				@Override
				public boolean test(SAMRecord rec) {
					if(rec.getReadUnmappedFlag()) return false;
					if(rec.isSecondaryOrSupplementary()) return false;
					if(rec.getDuplicateReadFlag()) return false;
					if(rec.getReadFailsVendorQualityCheckFlag()) return false;
					return true;
				}
			};
			
			for(SamReaderList samReader : samReaders) {
				SAMRecordIterator iter = samReader.get(0).iterator();
				while(iter.hasNext()) {
					final SAMRecord rec= iter.next();
					if(!samRecordFilter.test(rec)) continue;
					boolean skip=true;
					if(rec.getCigar()!=null && rec.getCigar().isClipped())
						{	
						skip=false;
						}
					
					if(skip) continue;
					
					if(	rec.getReadPairedFlag()) {
						if(!rec.getMateUnmappedFlag()) {
							SAMRecordIterator iter2 = samReader.get(1).query(
									rec.getMateReferenceName(),
									rec.getMateAlignmentStart(),
									rec.getMateAlignmentStart(),
									false
									);
							while(iter2.hasNext()) {
								final SAMRecord rec2 = iter2.next();
								if(!samRecordFilter.test(rec2)) continue;
								if(!rec2.getReadName().equals(rec.getReadName())) continue;
								if(rec2.getFirstOfPairFlag()==rec.getFirstOfPairFlag()) continue;
								if(rec2.getSecondOfPairFlag()==rec.getSecondOfPairFlag()) continue;
								challenge(rec,rec2);
								break;
								}
							iter2.close();
							}
						else
							{
							challenge(rec,null);
							}
						}
					}
				iter.close();
			}
			
			
			
			r =  null;
			out =  openFileOrStdoutAsPrintStream(outputFile);
			final SAMFileHeader header=r.getFileHeader();
			OtherCanonicalAlignFactory xpalignFactory=new OtherCanonicalAlignFactory(header);
			int prev_tid=-1;
			short coverage[]=null;
			short max_coverage=0;
			//w=swf.make(header, System.out);
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			final SAMRecordIterator it= r.iterator();
			for(;;)
				{
				SAMRecord rec=null;
				if(it.hasNext())
					{
					rec=progress.watch(it.next());
					if(rec.getReadUnmappedFlag()) continue;
					if(rec.isSecondaryOrSupplementary()) continue;
					if(rec.getDuplicateReadFlag()) continue;
					if(rec.getReadFailsVendorQualityCheckFlag()) continue;
					}
				
				if(rec==null || prev_tid==-1 || (prev_tid!=-1 && prev_tid!=rec.getReferenceIndex()))
					{
					if(coverage!=null)
						{

						while(max_coverage>=Math.max(1,this.min_coverage))
							{
							LOG.info("Scanning "+header.getSequence(prev_tid).getSequenceName()+" for cov:"+max_coverage);

							int chromStart0=0;
							while(chromStart0 < coverage.length)
								{
								if(coverage[chromStart0]==max_coverage)
									{
									coverage[chromStart0]=0;//reset that pos
									int chromEnd0=chromStart0+1;
									while(chromEnd0 < coverage.length && coverage[chromEnd0]==max_coverage)
										{
										coverage[chromEnd0]=0;//reset that pos
										++chromEnd0;
										}
									out.print(header.getSequence(prev_tid).getSequenceName());
									out.print('\t');
									out.print(chromStart0);
									out.print('\t');
									out.print(chromEnd0);
									out.print('\t');
									out.print(max_coverage);
									out.println();
									
									//reset 3'
									for(int x=chromEnd0;x<coverage.length && coverage[x]>0;++x)
										{
										coverage[x]=0;
										}
									//reset 5'
									for(int x=chromStart0-1;x>=0 && coverage[x]>0;--x)
										{
										coverage[x]=0;
										}
									chromStart0=chromEnd0;
									}
								else
									{
									++chromStart0;
									}
								}
							max_coverage--;
							}
						coverage=null;
						}
					
					if(rec==null) break;
					prev_tid=rec.getReferenceIndex();
					LOG.info("Alloc sizeof(short)*"+header.getSequence(prev_tid).getSequenceLength());
					coverage=new short[header.getSequence(prev_tid).getSequenceLength()];
					Arrays.fill(coverage,(short)0);
					max_coverage=0;
					}
				
				List<OtherCanonicalAlign> saList=xpalignFactory.getXPAligns(rec);
				if(saList.isEmpty()) continue;
				for(OtherCanonicalAlign xp:saList)
					{
					if(!xp.getReferenceName().equals(rec.getReferenceName())) continue;
					
					if(!rec.getReadNegativeStrandFlag()) //read is plus
						{
						if(!xp.getReadNegativeStrandFlag())
							{
							//info(xp+" vs strand "+rec.getReadNegativeStrandFlag());
							continue;//ignore both same strand
							}
						}
					else //read.strand=='-'
						{
						if(xp.getReadNegativeStrandFlag())
							{
							//info(xp+" vs strand "+rec.getReadNegativeStrandFlag());
							continue;//ignore both same strand
							}
						}
					if(Math.abs(rec.getUnclippedStart() - xp.getAlignmentStart()) > max_size_inversion) 
						{
						//info(xp+" vs pos "+rec.getUnclippedStart());
						continue;
						}
					int chromStart=Math.min(rec.getUnclippedStart(), xp.getAlignmentStart())-1;
					int chromEnd=Math.max(rec.getUnclippedEnd(), xp.getAlignmentStart())-1;
					for(int x=chromStart;x<=chromEnd && x< coverage.length;++x )
						{
						if(coverage[x]<Short.MAX_VALUE) coverage[x]++;
						if(max_coverage< coverage[x])
							{
							LOG.info("Max coverage "+max_coverage);
							max_coverage=coverage[x];
							}
						}
					}
				}
			
			it.close();
			r.close();r=null;
			progress.finish();
			return RETURN_OK;
			}

		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(r);
			CloserUtil.close(out);
			}
		}
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new SamShortInvertion().instanceMainWithExit(args);
		}
	}
