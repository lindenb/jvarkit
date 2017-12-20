package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.align.LongestCommonSequence;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

@Program(name="localrealignreads",description="Local Realignment of Reads")
public class LocalRealignReads extends Launcher
	{
	private static final Logger LOG = Logger.build(LocalRealignReads.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-R","--reference"},description="Indexed fasta Reference")
	private File faidxFile = null;
	
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs();
	

	private int MIN_ALIGN_LEN=15;
	
	
	
	
	private void align(
			final List<LongestCommonSequence.Hit> hits,
			final LongestCommonSequence matrix,
			final CharSequence S1,
		    int start1,
		    int end1,
			final CharSequence S2,
			int start2,
			int end2,
			int side
			)
		{
		if(end1-start1<MIN_ALIGN_LEN) return;
		if(end2-start2<MIN_ALIGN_LEN) return;
		
		LongestCommonSequence.Hit hit = matrix.align(
				S1,
				start1,end1,
				S2,
				start2,end2
				);
		int align_size = hit.size();
		if(align_size<MIN_ALIGN_LEN) return;
		hits.add(hit);
		
		if(side==-1 || side==1) align(hits,matrix,S1,start1,end1,S2,hit.getEndY(),end2,1);
		if(side==-1 || side==2) align(hits,matrix,S1,start1,end1,S2,0,hit.getStartY(),2);
		}
	
	@Override
	public int doWork(List<String> args) {
		if(this.faidxFile==null)
			{
			LOG.error("REFerence file missing;");
			return -1;
			}
		IndexedFastaSequenceFile  indexedFastaSequenceFile =null;
		SamReader samReader =null;
		SAMFileWriter w = null;
		SAMRecordIterator iter = null;
		GenomicSequence genomicSequence = null;
		LongestCommonSequence matrix =new LongestCommonSequence();
		try {
			indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidxFile);
			
			samReader  = openSamReader(oneFileOrNull(args));
			final SAMFileHeader header1 = samReader.getFileHeader();
			final SAMFileHeader header2 = header1.clone();
			header2.setSortOrder(SortOrder.unsorted);
			w = this.writingBamArgs.setReferenceFile(faidxFile).openSAMFileWriter(outputFile,header2, true);
			final SAMSequenceDictionaryProgress progress = new  SAMSequenceDictionaryProgress(header1);
			iter = samReader.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec = progress.watch(iter.next());
				if( rec.getReadUnmappedFlag() ||
					rec.isSecondaryOrSupplementary() ||
					rec.getReadFailsVendorQualityCheckFlag() ||
					rec.getDuplicateReadFlag()) {
					w.addAlignment(rec);
					continue;
					}
				
				final Cigar cigar = rec.getCigar(); 
				final int nCigarElement = cigar.numCigarElements();
				if(	nCigarElement <2) {
					w.addAlignment(rec);
					continue;
					}
				final  CigarElement ce5 =  cigar.getCigarElement(0);
				final  CigarElement ce3 =  cigar.getCigarElement(nCigarElement-1);
				if(!(
					(ce3.getOperator().equals(CigarOperator.S) && ce3.getLength()>=MIN_ALIGN_LEN) ||
					(ce5.getOperator().equals(CigarOperator.S) && ce5.getLength()>=MIN_ALIGN_LEN)) )
					{
					w.addAlignment(rec);
					continue;
					}
				
				if( genomicSequence == null || 
					!genomicSequence.getChrom().equals(rec.getReferenceName())
					)
					{
					genomicSequence = new GenomicSequence(indexedFastaSequenceFile, rec.getReferenceName());
					}
				
				final CharSequence readseq = rec.getReadString();
				
				
				/** short invert */
				if(ce5.getOperator()==CigarOperator.S && ce5.getLength()>=MIN_ALIGN_LEN)
					{
					CharSequence clipseq = new SubSequence(readseq, 0, ce5.getLength());
					CharSequence revcomp = new RevCompCharSequence(clipseq);
					LongestCommonSequence.Hit hit = matrix.align(genomicSequence,
							Math.max(0,rec.getUnclippedStart()-rec.getReadLength()),
							Math.min(rec.getUnclippedEnd()+rec.getReadLength(),genomicSequence.length()),
							revcomp,
							0,
							revcomp.length()
							);
					if(hit.size()>=MIN_ALIGN_LEN)
						{
						System.err.println("REVCOMP5' "+hit.getMatchingSequence()+" "+clipseq+" for "+readseq+" "+rec.getReadName()+" "+rec.getCigarString());
						}
					/*
					
					hit = matrix.align(
							readseq, 0, readseq.length(),
							revcomp,
							0,
							revcomp.length()
							);
					if(hit.size()>=MIN_ALIGN_LEN)
						{
						System.err.println("REVCOMP5' vs ITSELF: "+hit.getMatchingSequence()+" "+clipseq+" for "+readseq+" "+rec.getReadName()+" "+rec.getCigarString());
						}
					*/
					}
				if(ce3.getOperator()==CigarOperator.S && ce3.getLength()>=MIN_ALIGN_LEN)
					{
					CharSequence clipseq = new SubSequence(readseq,readseq.length()-ce3.getLength(),readseq.length());
					CharSequence revcomp = new RevCompCharSequence(clipseq);
					LongestCommonSequence.Hit hit = matrix.align(genomicSequence,
							Math.max(0,rec.getUnclippedStart()-rec.getReadLength()),
							Math.min(rec.getUnclippedEnd()+rec.getReadLength(),genomicSequence.length()),
							revcomp,
							0,
							revcomp.length()
							);
					if(hit.size()>=MIN_ALIGN_LEN)
						{
						System.err.println("REVCOMP3' "+hit.getMatchingSequence()+" "+clipseq+" for "+readseq+" "+rec.getReadName());
						}
					/*
					hit = matrix.align(
							readseq, 0, readseq.length(),
							revcomp,
							0,
							revcomp.length()
							);
					if(hit.size()>=MIN_ALIGN_LEN)
						{
						System.err.println("REVCOMP3' vs ITSELF: "+hit.getMatchingSequence()+" "+clipseq+" for "+readseq+" "+rec.getReadName()+" "+rec.getCigarString());
						}*/
					}
				
				/* other */
				List<LongestCommonSequence.Hit> hits = new ArrayList<>();
				align(hits, matrix,
						genomicSequence,
						Math.max(0,rec.getUnclippedStart()-rec.getReadLength()),
						Math.min(rec.getUnclippedEnd()+rec.getReadLength(),genomicSequence.length()),
						readseq,
						0,readseq.length(),
						-1
						);
				if(hits.size()<2000) continue;
				for(LongestCommonSequence.Hit hit:hits)
					{
					int readstart=0;
					boolean in_M=false;
					for(int i=0;i< nCigarElement;++i)
						{
						final CigarElement ce = cigar.getCigarElement(i);
						if(ce.getOperator().consumesReadBases())
							{
							int readend = readstart + ce.getLength(); 
							if( ce.getOperator() == CigarOperator.M ||
								ce.getOperator() == CigarOperator.X ||
								ce.getOperator() == CigarOperator.EQ)
								{
								if(!(hit.getEndY() <=readstart || readend<= hit.getStartY()))
									{
									in_M=true;
									break;
									}
								}
							readstart=readend;
							}
						}
					if(in_M) continue;
					
					int align_size =  hit.size();
					System.err.println(rec.getReadName()+" "+rec.getCigarString()+" "+rec.getAlignmentStart()+"-"+rec.getAlignmentEnd());
					System.err.println("REF: "+hit.getStartX()+"-"+hit.getEndX());
					System.err.println("READ "+hit.getStartY()+"-"+hit.getEndY());

					System.err.print("REF  :");
					for(int i=0;i< align_size;++i)
						{
						System.err.print(genomicSequence.charAt(hit.getStartX()+i));
						}
					System.err.println();
					System.err.print("READ :");
					
					for(int i=0;i< align_size;++i)
						{
						System.err.print(readseq.charAt(hit.getStartY()+i));
						}
					System.err.println();
					
					
					
						
						
					
					}
				
				System.err.println();
				
				}
			progress.finish();
			return RETURN_OK;
			}
		catch (Exception e) {
			return wrapException(e);
			}
		finally
			{
			genomicSequence=null;
			CloserUtil.close(iter);
			CloserUtil.close(samReader);
			CloserUtil.close(w);
			CloserUtil.close(indexedFastaSequenceFile);
			matrix=null;
			}
		}
	public static void main(String[] args) {
		new LocalRealignReads().instanceMainWithExit(args);
	}

}
