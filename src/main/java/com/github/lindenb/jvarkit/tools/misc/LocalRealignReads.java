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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

public class LocalRealignReads extends AbstractLocalRealignReads
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(LocalRealignReads.class);
	private int MIN_ALIGN_LEN=15;
	private static class Hit
		{
		int s1_from;
		int s1_to;
		int s2_from;
		int s2_to;
		}
	
	private static class Matrix
		{
		int width;
		int height;
		private int[] array=null;
	
		public Hit align(
				    final CharSequence S1,
				    int start1,
				    int end1,
					final CharSequence S2,
					int start2,
					int end2
					)
			{
			final int L1=(end1-start1);
			final int L2=(end2-start2);
			/** resize matrix */
			this.width = L1+1;
			this.height = L2+1;
			if( array==null || array.length< (width*height) )
				{
				this.array = new int[width*height];
				}
			/** reset matrix */
			Arrays.fill(this.array, 0);
			
			int best_x=0;
			int best_y=0;
			int max_score=0;
			this.array[0]=0;
			for(int x=0;x< L1 ;++x)
				{
				this.array[x+1]=0;
				}
			for(int y=0;y< L2 ;++y)
				{
				this.array[(y+1)*width]=0;
				for(int x=0;x< L1 ;++x)
					{
					char c1 = S1.charAt(start1 + x);
					char c2 = S2.charAt(start2 + y);
					int v;
					if( c1 == c2)
						{
						v = 1 + 
							this.array[(y)*width+(x)]//diagonal
							;
						}
					else
						{
						v = 0;
						}
					this.array[(y+1)*width+(x+1)] = v;
					if(v>max_score)
						{
						best_x  = x;
						best_y  = y;
						max_score=v;
						}
					}

				}
			Hit hit= new Hit();
			hit.s1_from = start1 + best_x - (max_score-1);
			hit.s1_to =   start1 + best_x+1;
			hit.s2_from = start2 + best_y - (max_score-1);
			hit.s2_to =   start2 + best_y+1;
			return hit;
			}
		}
	
	private void align(
			final List<Hit> hits,
			final Matrix matrix,
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
		
		Hit hit = matrix.align(
				S1,
				start1,end1,
				S2,
				start2,end2
				);
		int align_size = (hit.s2_to-hit.s2_from);
		if(align_size<MIN_ALIGN_LEN) return;
		hits.add(hit);
		
		if(side==-1 || side==1) align(hits,matrix,S1,start1,end1,S2,hit.s2_to,end2,1);
		if(side==-1 || side==2) align(hits,matrix,S1,start1,end1,S2,0,hit.s2_from,2);
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.faidxFile==null)
			{
			return wrapException("REFerence file missing;");
			}
		IndexedFastaSequenceFile  indexedFastaSequenceFile =null;
		SamReader samReader =null;
		SAMFileWriter w = null;
		SAMRecordIterator iter = null;
		GenomicSequence genomicSequence = null;
		Matrix matrix =new Matrix();
		try {
			indexedFastaSequenceFile = new IndexedFastaSequenceFile(this.faidxFile);
			
			samReader  = openSamReader(inputName);
			final SAMFileHeader header1 = samReader.getFileHeader();
			final SAMFileHeader header2 = header1.clone();
			header2.setSortOrder(SortOrder.unsorted);
			w = openSAMFileWriter(header2, true);
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
				
				List<Hit> hits = new ArrayList<>();
				align(hits, matrix,
						genomicSequence,
						Math.max(0,rec.getUnclippedStart()-rec.getReadLength()),
						Math.min(rec.getUnclippedEnd()+rec.getReadLength(),genomicSequence.length()),
						readseq,
						0,readseq.length(),
						-1
						);
				if(hits.size()<2) continue;
				for(Hit hit:hits)
					{
					int readstart=0;
					boolean in_M=false;
					for(int i=0;i< nCigarElement;++i)
						{
						final CigarElement ce = cigar.getCigarElement(i);
						if(ce.getOperator().consumesReadBases())
							{
							int readend = readstart + ce.getLength(); 
							if(ce.getOperator() == CigarOperator.M || ce.getOperator() == CigarOperator.X || ce.getOperator() == CigarOperator.EQ)
								{
								if(!(hit.s2_to <=readstart || readend<= hit.s2_from))
									{
									in_M=true;
									break;
									}
								}
							readstart=readend;
							}
						}
					if(in_M) continue;
					
					int align_size = (hit.s2_to-hit.s2_from);
					System.err.println(rec.getReadName()+" "+rec.getCigarString()+" "+rec.getAlignmentStart()+"-"+rec.getAlignmentEnd());
					System.err.println("REF: "+hit.s1_from+"-"+hit.s1_to);
					System.err.println("READ "+hit.s2_from+"-"+hit.s2_to);

					System.err.print("REF  :");
					for(int i=0;i< align_size;++i)
						{
						System.err.print(genomicSequence.charAt(hit.s1_from+i));
						}
					System.err.println();
					System.err.print("READ :");
					
					for(int i=0;i< align_size;++i)
						{
						System.err.print(readseq.charAt(hit.s2_from+i));
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
