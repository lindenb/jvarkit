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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.util.ucsc;

import java.io.BufferedReader;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.tribble.Feature;
import htsjdk.tribble.annotation.Strand;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.AbstractCharSequence;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.DelegateCharSequence;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;


public class KnownGene implements Iterable<Integer>,Feature
	{
	public static final String OPT_KNOWNGENE_DESC=
			"UCSC knownGene File/URL."
			+ " The knowGene format is a compact alternative to GFF/GTF because one transcript is described using only one line."
			+ "	Beware chromosome names are formatted the same as your REFERENCE. A typical KnownGene file is http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz ."
			+ "If you only have a gff file, you can try to generate a knownGene file with [http://lindenb.github.io/jvarkit/Gff2KnownGene.html](http://lindenb.github.io/jvarkit/Gff2KnownGene.html)";
	
	/** returns the UCSC URL for knownGene for the given UCSC build e.g: 'hg19' */
	public static String getUri(final String ucscBuild)
		{
		return "http://hgdownload.cse.ucsc.edu/goldenPath/"+ ucscBuild +"/database/wgEncodeGencodeBasicV19.txt.gz";
		}
	
	/** returns the default UCSC URL for knownGene */
	public static String getDefaultUri()
		{
		return getUri("hg19");
		}
	
	public enum PosType {
		INTERGENIC,
		UTR5,
		EXON,
		INTRON,
		UTR3
	}
	
	private String name;
	private String chrom;
	private char strand;
	private int txStart;
	private int txEnd;
	private int cdsStart;
	private int cdsEnd;
	private int exonStarts[];
	private int exonEnds[];
	private Map<String,String> attributes;
	
	
	@Override
	@Deprecated
	public final  String getChr() {
		return getContig();
		}
	
	@Override
	public final String getContig() {
		return getChromosome();
		}
	
	/** +1 based position. This is txStart+1 */
	@Override
	public final int getStart() {
		return getTxStart() + 1;
		}
	
	/** +1 based position. This is getTxEnd */
	@Override
	public final  int getEnd() {
		return getTxEnd();
		}
	
	/** returns true if cdsStart==cdsEnd */
	public boolean isNonCoding()
		{
		return getCdsStart()==getCdsEnd();
		}
	
	/** returns position type for given genomic position , never null*/
	public PosType getPositionTypeAt(int genomic0)
		{
		if(genomic0<this.getTxStart() || genomic0>=this.getTxEnd()) return PosType.INTERGENIC;
		if(genomic0 < this.getCdsStart()) return isPositiveStrand()?PosType.UTR5:PosType.UTR3;
		if(genomic0 >= this.getCdsEnd()) return isPositiveStrand()?PosType.UTR3:PosType.UTR5;
		for(int eidx=0;eidx< this.getExonCount();++eidx)
			{
			if(getExonStart(eidx)<=genomic0 && genomic0<getExonEnd(eidx)) return PosType.EXON;
			}
		return PosType.INTRON;
		}
	
	
	/** get transcript length (cumulative exons sizes )*/
	public int getTranscriptLength() {
		return getExons().stream().
				mapToInt(T->(T.getEnd()-T.getStart())).
				sum();
		}
	
	/** get UTR5 + UTR3 size */
	public int getUTRLength() {
		return getUTR5Length()+getUTR3Length();
	}
	
	/** get UTR5 size */
	public int getUTR5Length() {
		int n=0;
		for(int i=0;i< getExonCount();++i) {
			if(getExonEnd(i)< getCdsStart()) {
				n += getExonEnd(i)-getExonStart(i);
				}
			else if(getCdsStart()< getExonStart(i))
				{
				break;
				}
			else
				{
				n+= getCdsStart()-getExonStart(i);
				break;
				}
			}
		return n;
		}

	/** get UTR3 size */
	public int getUTR3Length() {
		int n=0;
		for(int i= getExonCount()-1;i>=0;--i) {
			if(getCdsEnd()<getExonStart(i) ) {
				n += getExonEnd(i)-getExonStart(i);
				}
			else if(getExonEnd(i)< getCdsEnd())
				{
				break;
				}
			else
				{
				n+= getExonEnd(i)-getCdsEnd();
				break;
				}
			}
		return n;
		}

	
	public abstract class AbstractSegment
		{
		/** shortcut to getGene().getContig() */
		public String getContig() {
			return this.getGene().getContig();
			}
		
		public KnownGene getGene()
			{
			return KnownGene.this;
			}
		
		public boolean isPositiveStrand()
	    	{
	    	return getGene().isPositiveStrand();
	    	}
	
		public boolean isNegativeStrand()
	    	{
	    	return getGene().isNegativeStrand();
	    	}
		public abstract String getName();
		/** the zero based position */
		public abstract int getStart();
		/** the zero based position */
		public abstract int getEnd();
		
		/** return true if the segment contans the genomic position 0-based */
		public boolean contains(int position0)
			{
			return getStart()<=position0 && position0< getEnd();
			}
		}
	
	public abstract class Segment 
		extends AbstractSegment 
		implements Iterable<Integer>
		{
		private final int index;
		protected Segment(int index)
			{
			this.index=index;
			}
		
		/** return 0 based index */
		public final int getIndex()
			{
			return index;
			}
		
		
		@Override
		public Iterator<Integer> iterator()
			{
			return iterator(false);
			}
		
		/** returns an iterator over all the 0-based genomic position of the KnownGene
		 * if useTranscriptDirection==true and strand is '-', will go from 3' to 5' (decreasing numbers)
		 * */
		public Iterator<Integer> iterator(boolean useTranscriptDirection)
			{
			final IntIter iter=new IntIter();
			if(useTranscriptDirection && this.isNegativeStrand())
				{
				iter.beg=this.getEnd()-1;
				iter.end=this.getStart()-1;
				iter.shift=-1;
				}
			else
				{
				iter.beg=this.getStart();
				iter.end=this.getEnd();
				iter.shift=1;
				}
			return iter;
			}
		
		
		public abstract boolean isSplicingAcceptor(int position);
		public abstract boolean isSplicingDonor(int position);
		public boolean isSplicing(int position)
			{
			return isSplicingAcceptor(position) || isSplicingDonor(position);
			}
		
		}
		
	public class Exon extends Segment
		{
		private Exon(final int index)
			{
			super(index);
			}
		/** return true if this exon does NOT overlap a CDS = exon in UTR or gene is non-coding*/
		public boolean isNonCoding()
			{
			if(getGene().isNonCoding()) return true;
			if(this.getEnd()<=getGene().getCdsStart()) return true;
			if(getGene().getCdsEnd()<=this.getStart()) return true;
			return false;
			}
		@Override
		public String getName()
			{
			if(getGene().isPositiveStrand())
				{
				return "Exon "+(getIndex()+1);
				}
			else
				{
				return "Exon "+(getGene().getExonCount()-getIndex());
				}
			}
		
		@Override
		public int getStart()
			{
			return getGene().getExonStart(getIndex());
			}
		
		@Override
		public int getEnd()
			{
			return getGene().getExonEnd(getIndex());
			}
		
		@Override
		public String toString()
			{
			return getName();
			}
		
		/** return null if there is no next intron */
		public Intron getNextIntron()
			{
			if(getIndex()+1>=getGene().getExonCount()) return null;
			return getGene().getIntron(getIndex());
			}
		/** return null if there is no prev intron */
		public Intron getPrevIntron()
			{
			if(getIndex()<=0) return null;
			return getGene().getIntron(getIndex()-1);
			}
		
		@Override
		public boolean isSplicingAcceptor(int position)
			{
			if(!contains(position)) return false;
			if(isPositiveStrand())
				{
				if(getIndex()== 0) return false;
				return position==getStart();
				}
			else
				{
				if(getIndex()+1== getGene().getExonCount()) return false;
				return position==getEnd()-1;
				}
			}
		
		@Override
		public boolean isSplicingDonor(int position)
			{
			if(!contains(position)) return false;
			if(isPositiveStrand())
				{
				if(getIndex()+1== getGene().getExonCount()) return false;
				return  (position==getEnd()-1) ||
						(position==getEnd()-2) ||
						(position==getEnd()-3) ;
				}
			else
				{
				if(getIndex()== 0) return false;
				return  (position==getStart()+0) ||
						(position==getStart()+1) ||
						(position==getStart()+2) ;
				}
			}
		/** return true if it's the first exon in the transcription direction */
		public boolean isFirstExon()
			{
			if(this.getGene().isPositiveStrand() &&  this.getIndex()==0) return true;
			if(this.getGene().isNegativeStrand() && this.getIndex()+1 == this.getGene().getExonCount()) return true;
			return false;
			}
		/** return true if it's the last exon in the transcription direction */
		public boolean isLastExon()
			{
			if(this.getGene().isNegativeStrand() &&  this.getIndex()==0) return true;
			if(this.getGene().isPositiveStrand() && this.getIndex()+1 == this.getGene().getExonCount()) return true;
			return false;
			}
		}
		
	/** an intron in a knownGene */
	public class Intron extends Segment
			{
			Intron(int index)
				{
				super(index);
				}
			/** the zero based position */
			@Override
			public int getStart()
				{
				return getGene().getIntronStart(getIndex());
				}
			/** the zero based position */
			@Override
			public int getEnd()
				{
				return getGene().getIntronEnd(getIndex());
				}
			
			@Override
			/** get intron name using one-based human order (according to strand */
			public String getName() {
				if(getGene().isPositiveStrand())
					{
					return "Intron "+(getIndex()+1);
					}
				else
					{
					return "Intron "+(getGene().getExonCount()-getIndex());
					}
				}

			public boolean isSplicingAcceptor(int position)
				{
				if(!contains(position)) return false;
				if(isPositiveStrand())
					{
					return  (position==getEnd()-1) ||
							(position==getEnd()-2);
					}
				else
					{
					return	position==getStart() ||
							position==getStart()+1;
					}
				}
			

			public boolean isSplicingDonor(int position)
				{
				if(!contains(position)) return false;
				if(isPositiveStrand())
					{
					return	position==getStart() ||
							position==getStart()+1;
							
					}
				else
					{
					return  (position==getEnd()-1) ||
							(position==getEnd()-2);
					}
				}
			
			}
	
		/**
		 * 
		 * KnownGene 
		 * 
		 */
		public KnownGene()
			{
			this.name="";
			this.chrom="";
			this.strand='+';
			this.txStart=0;
			this.txEnd=0;
			this.cdsStart=0;
			this.cdsEnd=0;
			this.exonStarts=new int[]{0};
			this.exonEnds=new int[]{0};
			}
		/** column to be retrieved in sql query */
		public static final String SQL_COLUMNS[]={
			"name","chrom","strand",
			"txStart","txEnd",
			"cdsStart","cdsEnd",
			"exonCount",
			"exonStarts",
			"exonEnds"
			};
		
		public KnownGene(final ResultSet row) throws SQLException
			{
			this.name = row.getString( SQL_COLUMNS[0]);
			this.chrom= row.getString(SQL_COLUMNS[1]);
	        this.strand = row.getString(SQL_COLUMNS[2]).charAt(0);
	        this.txStart = row.getInt(SQL_COLUMNS[3]);
	        this.txEnd = row.getInt(SQL_COLUMNS[4]);
	        this.cdsStart= row.getInt(SQL_COLUMNS[5]);
	        this.cdsEnd= row.getInt(SQL_COLUMNS[6]);
	        int exonCount=row.getInt(SQL_COLUMNS[7]);
	        this.exonStarts = new int[exonCount];
	        this.exonEnds = new int[exonCount];
	        
            int index=0;
            for(String s:row.getString(SQL_COLUMNS[8]).split("[,]"))
            	{
            	this.exonStarts[index++]=Integer.parseInt(s);
            	}
            index=0;
            for(String s:row.getString(SQL_COLUMNS[9]).split("[,]"))
            	{
            	this.exonEnds[index++]=Integer.parseInt(s);
            	}

			}
		
		private static final CharSplitter COMMA = CharSplitter.COMMA;
		
		public KnownGene(final String tokens[])
			{
			/* first column may be the bin column 
			 * we use the strand to detect the format */
			final int binIdx=tokens[2].equals("+") || tokens[2].equals("-")?0:1;
			
			this.name = tokens[binIdx + 0];
			this.chrom= tokens[binIdx + 1];
	        this.strand = tokens[binIdx + 2].charAt(0);
	        this.txStart = Integer.parseInt(tokens[binIdx + 3]);
	        this.txEnd = Integer.parseInt(tokens[binIdx + 4]);
	        this.cdsStart= Integer.parseInt(tokens[binIdx + 5]);
	        this.cdsEnd= Integer.parseInt(tokens[binIdx + 6]);
	        final int exonCount=Integer.parseInt(tokens[binIdx + 7]);
	        this.exonStarts = new int[exonCount];
	        this.exonEnds = new int[exonCount];
	        
            int index=0;
            for(final String s: COMMA.splitAsStringList(tokens[binIdx + 8]))
            	{
            	this.exonStarts[index++]=Integer.parseInt(s);
            	}
            index=0;
            for(final String s: COMMA.splitAsStringList(tokens[binIdx + 9]))
            	{
            	this.exonEnds[index++]=Integer.parseInt(s);
            	}
			}
		
		/** returns knownGene ID */
		public String getName()
			{
			return this.name;
			}
		
		public void setName(final String name) {
			this.name = name;
		}
		
		/** returns chromosome name */
		public String getChromosome()
			{
			return this.chrom;
			}
		
		public void setChrom(final String chrom) {
			this.chrom = chrom;
			}
		
		/** returns the strand */
		public Strand getStrand()
			{
			switch(strand)
				{
				case '+': return Strand.POSITIVE;
				case '-': return Strand.NEGATIVE;
				default: return Strand.NONE;
				}
			}
		
		public void setStrand(char strand)
			{
			this.strand = strand;
			}
		
		public boolean isPositiveStrand()
        	{
        	return getStrand()==Strand.POSITIVE;
        	}

		public boolean isNegativeStrand()
	    	{
	    	return getStrand()==Strand.NEGATIVE;
	    	}

		public void setTxStart(final int txStart) {
			this.txStart = txStart;
		}
		public int getTxStart()
			{
			return this.txStart;
			}
		
		public void setTxEnd(final int txEnd) {
			this.txEnd = txEnd;
			}

		public int getTxEnd()
			{
			return this.txEnd;
			}
		
		public void setCdsStart(final int cdsStart) {
			this.cdsStart = cdsStart;
			}
		public int getCdsStart()
			{
			return this.cdsStart;
			}
		
		public void setCdsEnd(final int cdsEnd) {
			this.cdsEnd = cdsEnd;
			}
		public int getCdsEnd()
			{
			return this.cdsEnd;
			}
		

		public void setExonBounds(final int exonCount,final String exonStarts,final String exonEnds)
			{
			this.exonStarts=new int[exonCount];
			this.exonEnds=new int[exonCount];
			int i=0;
			for(final String s: COMMA.split(exonStarts))
				{
				this.exonStarts[i++]=Integer.parseInt(s);
				}
			i=0;
			for(final String s: COMMA.split(exonEnds))
				{
				this.exonEnds[i++]=Integer.parseInt(s);
				}
			}
		
		public int getExonStart(int index)
			{
			return this.exonStarts[index];
			}
		

		public int getExonEnd(int index)
			{
			return this.exonEnds[index];
			}
		
		/** get 0-based start for the index-th intron */
		public int getIntronStart(int index)
			{
			return this.getExonEnd(index);
			}
		
		/** get 0-based end for the index-th intron */
		public int getIntronEnd(int index)
			{
			return this.getExonStart(index+1);
			}

		
		/** exon[0] is first from 5' to 3' on genomic */
		public Exon getExon(int index)
			{
			return new Exon(index);
			}
		/** intron[0] is intron after exon[0] from 5' to 3' on genomic */
		public Intron getIntron(int i)
			{
			return new Intron(i);
			}
		/** get number of exons */
		public int getExonCount()
			{
			return this.exonStarts.length;
			}
		/** get number of introns */
		public int getIntronCount()
			{
			return Math.max(getExonCount()-1,0);
			}
		
		public Map<String,String> getAttributes()
			{
			if(this.attributes==null) this.attributes=new HashMap<String, String>();
			return this.attributes;
			}
		
		
		public List<Exon> getExons()
			{
			final List<Exon> L=new ArrayList<>(getExonCount());
			for(int i=0;i< getExonCount();++i)
				{	
				L.add(getExon(i));
				}
			return L;
			}
		/** return a list of intron for this gene */
		public List<Intron> getIntrons()
			{
			final List<Intron> L=new ArrayList<>(getIntronCount());
			for(int i=0;i< getIntronCount();++i)
				{	
				L.add(getIntron(i));
				}
			return L;
			}
		
		@Override
		public Iterator<Integer> iterator()
			{
			return iterator(false);
			}
		
		/** returns an iterator over all the 0-based genomic position of the KnownGene
		 * if useTranscriptDirection==true and strand is '-', will go from 3' to 5' (decreasing numbers)
		 * */
		public Iterator<Integer> iterator(boolean useTranscriptDirection)
			{
			final IntIter iter=new IntIter();
			if(useTranscriptDirection && isNegativeStrand())
				{
				iter.beg=txEnd-1;
				iter.end=txStart-1;
				iter.shift=-1;
				}
			else
				{
				iter.beg=txStart;
				iter.end=txEnd;
				iter.shift=1;
				}
			return iter;
			}
		
		private static class IntIter implements Iterator<Integer>
			{
			int beg;
			int end;
			int shift;
			@Override
			public boolean hasNext() {
				return this.beg!=this.end;
				}
			@Override
			public Integer next()
				{
				int n=beg;
				beg+=shift;
				return n;
				}
			}
		
		
		
		
		/** an abstract RNA : start() and end() on genomic sequence need to be defined*/
		abstract class RNA extends DelegateCharSequence
			{
			/** cache length */
			private Integer _length=null;
			RNA(final CharSequence sequence)
				{
				super(sequence);
				}
			/** get Gene associated to this RNA */
			public final KnownGene getKnownGene()
				{
				return KnownGene.this;
				}
			/** start of mRNA (could be transcription or traduction */
			protected abstract int start();
			/** end of mRNA (could be transcription or traduction */
			protected abstract int end();
			
			public String getContig() {
				return getKnownGene().getContig();
				}
			
			/** return 0-based position in genomic coordinate !!! */
			public int getGenomicStart0() {
				return start();
			}
			/** return 0-based position in genomic coordinate !!! */
			public int getGenomicEnd0() {
				return end();
			}
			
		
			/** convert the genomic position to the position in the RNA, return -1 if RNA not in genomic pos */
			public int convertGenomicToRnaCoordinate(int genomicPos0)
				{
				int rnaPos0=0;
				if(getKnownGene().isPositiveStrand())
					{
					for(final Exon ex:getKnownGene().getExons())
						{
						if(this.start()>=ex.getEnd()) continue;
						if(this.end()<=ex.getStart()) break;
						final int beg=Math.max(this.start(), ex.getStart());
						final int end=Math.min(this.end(), ex.getEnd());
						if(beg<= genomicPos0 && genomicPos0<end)
							{
							return (genomicPos0-beg)+rnaPos0;
							}
						rnaPos0+=(end-beg);
						}
					}
				else
					{
					for(int idx=getKnownGene().getExonCount()-1;idx>=0;idx--)
						{
						final Exon ex=getExon(idx);
						if(this.start()>=ex.getEnd()) break;
						if(this.end()<=ex.getStart()) continue;
						final int beg=Math.max(this.start(), ex.getStart());
						final int end=Math.min(this.end(), ex.getEnd());
						if(beg<= genomicPos0 && genomicPos0<end)
							{
							return ((end-1)-genomicPos0) + rnaPos0;
							}
						rnaPos0+=(end-beg);
						}
					}
				return -1;
				}
			
			public int convertToGenomicCoordinate(int rnaPos0)
				{
				if(rnaPos0<0) throw new IllegalArgumentException("negative index:"+rnaPos0);
				if(rnaPos0>=this.length()) throw new IndexOutOfBoundsException("out of bound index:"+rnaPos0+"<"+this.length());
				
				if(getKnownGene().isPositiveStrand())
					{
					for(Exon ex:getKnownGene().getExons())
						{
						if(this.start()>=ex.getEnd()) continue;
						if(this.end()<=ex.getStart()) break;
						int beg=Math.max(this.start(), ex.getStart());
						int end=Math.min(this.end(), ex.getEnd());
						int len=end-beg;
						if(rnaPos0<len)
							{
							return beg+rnaPos0;
							}
						rnaPos0-=len;
						}
					}
				else
					{
					for(int idx=getKnownGene().getExonCount()-1;idx>=0;idx--)
						{
						Exon ex=getExon(idx);
						if(this.start()>=ex.getEnd()) break;
						if(this.end()<=ex.getStart()) continue;
						int beg=Math.max(this.start(), ex.getStart());
						int end=Math.min(this.end(), ex.getEnd());
						int len=end-beg;
						if(rnaPos0<len)
							{
							return (end-1)-rnaPos0;
							}
						rnaPos0-=len;
						}
					}
				return -1;
				}
			
			public Exon getExonAt(int rnaPos0)
				{
				for(final Exon ex:getKnownGene().getExons())
					{
					if(this.start()>=ex.getEnd()) continue;
					if(this.end()<=ex.getStart()) break;
					int beg=Math.max(this.start(), ex.getStart());
					int end=Math.min(this.end(), ex.getEnd());
					int n=(end-beg);
					if(rnaPos0<n) return ex;
					}
				return null;
				}
			
			@Override
			public int length()
				{
				if(this._length==null)
					{
					this._length=0;
					for(final Exon ex:getKnownGene().getExons())
						{
						if(this.start()>=ex.getEnd()) continue;
						if(this.end()<=ex.getStart()) break;
						int beg=Math.max(this.start(), ex.getStart());
						int end=Math.min(this.end(), ex.getEnd());
						this._length+=(end-beg);
						}
					}
				return this._length;
				}
			@Override
			public char charAt(int index0)
				{
				if(index0<0) throw new IllegalArgumentException("negative index:"+index0);
				if(index0>=this.length()) throw new IndexOutOfBoundsException("index:"+index0 +" < "+this.length());
				int n= convertToGenomicCoordinate(index0);
				if(n==-1) 	throw new IndexOutOfBoundsException("0<=index:="+index0+"<"+length());
				if(getKnownGene().isPositiveStrand())
					{
					return getDelegate().charAt(n);	
					}
				else
					{	
					return AcidNucleics.complement(getDelegate().charAt(n));	
					}
				}
			
			}
		
		
		/** same as RNA but without the UTR, starts is cdsStart, end is cdsEnd */
		public class CodingRNA  extends RNA
			{
			CodingRNA(final CharSequence sequence)
				{
				super(sequence);
				}
			@Override
			protected  int start()
				{
				return getKnownGene().getCdsStart();	
				}
			@Override
			protected  int end()
				{
				return getKnownGene().getCdsEnd();	
				}
			public Peptide getPeptide()
				{
				return new Peptide(GeneticCode.getStandard(), this);
				}
			}
		
		public CodingRNA getCodingRNA()
			{
			return getCodingRNA(new AbstractCharSequence()
				{
				@Override
				public int length() {
					return KnownGene.this.getTxEnd()+1;
				}
				
				@Override
				public char charAt(int index) {
					return 'N';
				}
				});
			}
		/** returns the coding RNA for this gene */
		public CodingRNA getCodingRNA(final CharSequence genomic)
			{
			return new CodingRNA(genomic);
			}
		
		
		public class MessengerRNA  extends RNA
			{
			MessengerRNA(final CharSequence sequence)
				{
				super(sequence);
				}
			@Override
			protected  int start()
				{
				return getKnownGene().getTxStart();	
				}
			@Override
			protected  int end()
				{
				return getKnownGene().getTxEnd();	
				}
			}
		
		public class Peptide extends DelegateCharSequence
			{
			private Integer _length=null;
			private final GeneticCode geneticCode;
			Peptide(final GeneticCode gc,CodingRNA rna)
				{
				super(rna);
				this.geneticCode=gc;
				}
			public CodingRNA getCodingRNA()
				{
				return CodingRNA.class.cast(getDelegate());
				}
			@Override
			public int length()
				{
				if(_length==null)
					{
					_length = getDelegate().length()/3;
					if(_length==0)
						{
						return 0;
						}
					final int idx1=getCodingRNA().convertToGenomicCoordinate((_length-1)*3+0);
					final int idx2=getCodingRNA().convertToGenomicCoordinate((_length-1)*3+1);
					final int idx3=getCodingRNA().convertToGenomicCoordinate((_length-1)*3+2);
					if(idx1==-1 || idx2==-1 || idx3==-1)
						{
						System.err.println("Bizarre pour "+KnownGene.this.getName()+" "+getStrand());
						return _length;
						}
					
					final char last=this.getGeneticCode().translate(
							getCodingRNA().charAt((_length-1)*3+0),	
							getCodingRNA().charAt((_length-1)*3+1),	
							getCodingRNA().charAt((_length-1)*3+2)
							);
					if(!Character.isLetter(last)) _length--;
					}
				return _length;
				}
			
			public char[] getCodon(int pepPos0)
				{
				return new char[]{
					getCodingRNA().charAt(pepPos0*3+0),	
					getCodingRNA().charAt(pepPos0*3+1),	
					getCodingRNA().charAt(pepPos0*3+2)
					};
				}

			
			/** convert the genomic position to the position in the peptide, return -1 if peptide not in genomic pos */
			public int convertGenomicToPeptideCoordinate(int genomicPos0)
				{
				int rnaIdx=getCodingRNA().convertGenomicToRnaCoordinate(genomicPos0);
				if(rnaIdx==-1) return -1;
				return rnaIdx/3;
				}
			
			
			
			public int[] convertToGenomicCoordinates(int pepPos0)
				{
				if(pepPos0<0) throw new IndexOutOfBoundsException("negative offset : "+pepPos0);
				if(pepPos0>=this.length()) throw new IndexOutOfBoundsException(" idx="+pepPos0+" and length="+this.length());
				
				return new int[]{
					getCodingRNA().convertToGenomicCoordinate(pepPos0*3+0),	
					getCodingRNA().convertToGenomicCoordinate(pepPos0*3+1),	
					getCodingRNA().convertToGenomicCoordinate(pepPos0*3+2)
					};
				}
			
			@Override
			public char charAt(final int pepPos0) {
				return this.getGeneticCode().translate(
						getCodingRNA().charAt(pepPos0*3+0),	
						getCodingRNA().charAt(pepPos0*3+1),	
						getCodingRNA().charAt(pepPos0*3+2)
						);
				}
			
			public GeneticCode getGeneticCode() {
				return this.geneticCode;
				}
			}
		
		@Override
		public String toString() {
			return getName()+"["+getContig()+":"+getTxStart()+"-"+getTxEnd()+"]";
			}
		
		
		/** load knownGene file/uri as an IntervalTreeMap. Intervals in the IntervalTreeMap are *1-based* (interval.start= kg.txStart+1)*/
		public static IntervalTreeMap<List<KnownGene>> loadUriAsIntervalTreeMap(
				final String uri,
				final Predicate<KnownGene> filterOrNull
				) throws IOException
			{
			final IntervalTreeMap<List<KnownGene>> treeMap = new IntervalTreeMap<>();
			BufferedReader in=null;
			try {
				in = IOUtils.openURIForBufferedReading(uri);
				String line;
				final CharSplitter tab = CharSplitter.TAB;
				while ((line = in.readLine()) != null) {
					if (line.isEmpty())
						continue;
					final String tokens[] = tab.split(line);
					final KnownGene g = new KnownGene(tokens);
					if(filterOrNull!=null && !filterOrNull.test(g)) continue;
					final Interval interval = new Interval(g.getContig(),g.getTxStart()+1,g.getTxEnd(),g.isNegativeStrand(),g.getName());
					List<KnownGene> L=  treeMap.get(interval);
					if(L==null) {
						L=new ArrayList<>(2);
						treeMap.put(interval, L);
						}
					
					L.add(g);
				}
				in.close();
				in = null;
				return treeMap;
				}
			finally {
				CloserUtil.close(in);
				}
			}
	}
