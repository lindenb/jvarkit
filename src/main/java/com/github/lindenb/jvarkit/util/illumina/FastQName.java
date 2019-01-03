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
package com.github.lindenb.jvarkit.util.illumina;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Optional;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.StringUtil;

import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;


/** FastQName */
public abstract class FastQName
	implements Comparable<FastQName>
	{
	
	private static final String SUFFIX=".fastq.gz";
	private static final int NO_LANE=-1;
	private static final int NO_SPLIT=-1;
	
	public enum Side {
		None,//first in list
		Forward,Reverse
		}
	protected Long nReads=null;
	protected FastQName()
		{
		}
	
	public abstract File getFile();
	public abstract boolean isUndetermined();
	public abstract String getSample();
	public boolean hasLane() { return getLane()!=NO_LANE;}
	public abstract int getLane();
	public boolean hasSplit() { return getSplit()!=NO_SPLIT;}
	public abstract int getSplit();
	public abstract Side getSide();
	public boolean hasSeqIndex() {
		return !StringUtil.isBlank(getSeqIndex());
	}
	public abstract String getSeqIndex();
	public abstract boolean isValid();
	
	

	
	private static class FastqNameBcl2Fastq2 extends FastQName {
		private File file;
		private int sample_index_in_samplesheet = -1;
		private int split = NO_SPLIT;
		private Side side = Side.None;
		private String sample;
		@Override
		public File getFile() {
			return file;
		}

		@Override
		public boolean isUndetermined() {
			return getSample().equals("Undetermined") ||
					sample_index_in_samplesheet == 0;
		}

		@Override
		public String getSample() {
			return sample;
		}

		@Override
		public int getLane() {
			return NO_LANE;
		}

		@Override
		public int getSplit() {
			return this.split;
		}

		@Override
		public Side getSide() {
			return side;
		}

		@Override
		public String getSeqIndex() {
			return null;
		}

		@Override
		public boolean isValid() {
			return true;
			}
		
		}
	
	/** parse as the new version of illumina bcl2fastq2 .eg: Pool2_S2_R2_001.fastq.gz*/
	public static Optional<FastQName> parseAsBcl2fastq2(final File file) {
		String fname = file.getName();
		if(!fname.endsWith(SUFFIX)) return Optional.empty();
		fname = fname.substring(0,fname.length()-SUFFIX.length());
		final String tokens[]=fname.split("[_]");
		if(tokens.length<4) return Optional.empty();
		
		final FastqNameBcl2Fastq2 fq = new FastqNameBcl2Fastq2();
		fq.file = file;
		/* parse 'split' */
		int column = tokens.length-1;
		try {
			fq.split = Integer.parseInt(tokens[column]);
		} catch(final NumberFormatException err) {
			return Optional.empty();
		}
		if(fq.split<0) return Optional.empty();
		/* parse 'side' */
		column--;
		if(tokens[column].equals("R1"))
			{
			fq.side = Side.Forward;
			}
		else if(tokens[column].equals("R2"))
			{
			fq.side = Side.Reverse;
			}
		else
			{
			return Optional.empty();
			}
		/* parse 'sample_index' */
		column--;
		if(!tokens[column].startsWith("S"))
			{
			return Optional.empty();
			}
		try {
			fq.sample_index_in_samplesheet = Integer.parseInt(tokens[column].substring(1));
		} catch(final NumberFormatException err) {
			return Optional.empty();
		}
		if(fq.sample_index_in_samplesheet<0) //0==undetermined
			{
			return Optional.empty();
			}
		
		fq.sample = String.join("_", Arrays.asList(tokens).subList(0, column));
		
		return Optional.of(fq);
	}
	
	private static class FastqNameImpl extends FastQName {
		protected File file;
		protected boolean is_valid=false;
		protected String seqIndex="";
		protected int lane = NO_LANE;
		protected Side side=Side.None;
		protected String sample=null;
		protected int split=NO_SPLIT;
		protected boolean undetermined=false;
	
		
		@Override
		public boolean isValid()
			{
			return is_valid;
			}
		
		@Override
		public File getFile() {
		return file;
		}
		
		@Override
		public boolean isUndetermined() {
		return undetermined;
		}
		
		@Override
		public String getSample() {
			return sample;
			}
		
		@Override
		public int getLane() {
			return lane;
			}
		
		@Override
		public int getSplit()
			{
			return split;
			}
		
		@Override
		public Side getSide()
			{
			return side;
			}
		
		@Override
		public String getSeqIndex()
			{
			return seqIndex;
			}
		
		}
	
	public static FastQName parse(final File f)
		{
		
		final Optional<FastQName> opt = parseAsBcl2fastq2(f);
		if(opt.isPresent()) return opt.get();
		
		FastqNameImpl fq=new FastqNameImpl();
		fq.file=f;
		
		if(!f.getName().endsWith(SUFFIX))
			{
			return fq;
			}
	
		String tokens[]=f.getName().split("[_]");
		//remove suffix
		int token_index=tokens.length-1;
		tokens[token_index]=tokens[token_index].substring(0,tokens[token_index].length()-SUFFIX.length());
		try
			{
			fq.split=Integer.parseInt(tokens[token_index]);
			}
		catch(Exception err)
			{
			return fq;
			}
		
		token_index--;
		if(token_index<0) return fq;
		

		if(tokens[token_index].equals("R1"))
			{
			fq.side=Side.Forward;
			}
		else if(tokens[token_index].equals("R2"))
			{
			fq.side=Side.Reverse;
			}
		else
			{
			return fq;
			}
		
		token_index--;
		if(token_index<0) return fq;
		
		
		if(!tokens[token_index].startsWith("L"))
			{
			return fq;
			}
		try
			{
			fq.lane=Integer.parseInt(tokens[token_index].substring(1));
			}
		catch(Exception err)
			{
			return fq;
			}
		
		token_index--;
		if(token_index<0) return fq;
		
		
		if(tokens[token_index].equalsIgnoreCase("undetermined"))
			{
			fq.is_valid=true;
			fq.undetermined=true;
			return fq;
			}
		/** update 2017-11-28, index can contain an hyphen for double indexing eg. SX1R_1_CTGAAGCT-AGGCGAAG_L001_R1_001.fastq.gz*/
		if(token_index>0 && 
				(
				 tokens[token_index].matches("[ATGC]{4,9}") ||
				 tokens[token_index].matches("[ATGC]{4,9}\\-[ATGC]{4,9}")
				 ))
			{
			fq.seqIndex=tokens[token_index];
			token_index--;
			}
		else if(token_index>0 && tokens[token_index].matches("S[0-9]+"))
			{
			//MiSeq
			token_index--;
			}
			
		
		final StringBuilder b=new StringBuilder();
		for(int i=0;i<=token_index;++i)
			{
			if(i>0) b.append('_');
			b.append(tokens[i]);
			}
		
		fq.sample=b.toString();
		fq.is_valid=true;
		return fq;
		}
	
	
	
	
	
	public long countReads() throws IOException
		{
		if(nReads==null)
			{
			long n=0L;
			FastqReader r=null;
			try
				{
				r=new FourLinesFastqReader(this.getFile());
				r.setValidationStringency(ValidationStringency.LENIENT);
				while(r.hasNext())
					{
					r.next();
					++n;
					}
				}
			catch(Exception err)
				{
				nReads=-1L;
				throw new IOException(err);
				}
			finally
				{
				if(r!=null)r.close();
				}
			
			nReads=n;
			}
		return nReads;
		}
	
	/** test wether the file contains at least one FASTQ record */
	public boolean isEmpty() throws IOException
		{
		int nLines=0;
		int c;
		GZIPInputStream in=null;
			try
				{
				in = new GZIPInputStream( new FileInputStream(getFile()));
				while((c=in.read())!=-1)
					{
					if(c=='\n')
						{
						nLines++;
						if(nLines>3) break;
						}
					}
				return nLines<4;
				}
		finally
			{
			if(in!=null) try{ in.close();}catch(IOException err2){}
			}
		}

	


	@Override
	public String toString() {
		return "FastQName [file=" + getFile() + ", is_valid=" + isValid()
				+ ", seqIndex=" + getSeqIndex() + ", lane=" + getLane() + ", side="
				+ getSide() + ", sample=" + getSample() + ", split=" + getSplit()
				+ ", undetermined=" + isUndetermined() + "]";
	}


	@Override
	public int hashCode() {
		return getFile().hashCode();
		}
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==null || !(obj instanceof FastQName)) return false;
		if(obj==this ) return true;
		return getFile().equals(((FastQName)obj).getFile());
		}
	
	@Override
	public int compareTo(final FastQName o) {
		return getFile().getPath().compareTo(o.getFile().getPath());
		}
	
	public boolean isComplementOf(final FastQName other)
		{
		if(!isValid() || !other.isValid()) return false;
		if(! ((this.isUndetermined() && other.isUndetermined()) || (this.getSample().equals(other.getSample()))) )
			{
			return false;
			}
		
		if(!this.isUndetermined() &&
				this.hasSeqIndex() &&
				other.hasSeqIndex() &&
				!this.getSeqIndex().equals(other.getSeqIndex()))
			{
			return false;
			}
		if(this.hasLane() && other.hasLane() && this.getLane()!=other.getLane()) return false;
		if(this.hasSplit() && other.hasSplit() && this.getSplit()!=other.getSplit()) return false;
		if(this.hasSeqIndex() && other.hasSeqIndex() &&  !this.getSeqIndex().equals(other.getSeqIndex())  ) return false;
		if(this.getSide()==other.getSide()) return false;
		return true;
		}
    		
	
	}
