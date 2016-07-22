/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.ValidationStringency;

import com.github.lindenb.jvarkit.util.picard.FastqReader;
import com.github.lindenb.jvarkit.util.picard.FourLinesFastqReader;


/** FastQName */
public class FastQName
	implements Comparable<FastQName>
	{
	private static final String SUFFIX=".fastq.gz";
	private File file;
	private boolean is_valid=false;
	private String seqIndex="";
	private int lane;
	private Side side=Side.None;
	private String sample=null;
	private int split=-1;
	private boolean undetermined=false;
	private Long nReads=null;
	
	public enum Side {
		None,//first in list
		Forward,Reverse
		}
	
	private FastQName()
		{
		}
	
	public boolean isUndetermined()
		{
		return undetermined;
		}
	
	public String getSample() {
		return sample;
		}
	
	public int getLane() {
		return lane;
		}
	
	public int getSplit()
		{
		return split;
		}
	
	public Side getSide()
		{
		return side;
		}
	
	public String getSeqIndex()
		{
		return seqIndex;
		}
	
	
	
	public static FastQName parse(File f)
		{
		
		FastQName fq=new FastQName();
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
		
		if(token_index>0 && tokens[token_index].matches("[ATGC]{4,9}"))
			{
			fq.seqIndex=tokens[token_index];
			token_index--;
			}
		else if(token_index>0 && tokens[token_index].matches("S[0-9]+"))
			{
			//MiSeq
			token_index--;
			}
			
		
		StringBuilder b=new StringBuilder();
		for(int i=0;i<=token_index;++i)
			{
			if(i>0) b.append('_');
			b.append(tokens[i]);
			}
		
		fq.sample=b.toString();
		fq.is_valid=true;
		return fq;
		}
	
	
	public File getFile()
		{
		return file;
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

	
	/** the record was parsed */
	public boolean isValid()
		{
		return is_valid;
		}

	@Override
	public String toString() {
		return "FastQName [file=" + file + ", is_valid=" + is_valid
				+ ", seqIndex=" + seqIndex + ", lane=" + lane + ", side="
				+ side + ", sample=" + sample + ", split=" + split
				+ ", undetermined=" + undetermined + "]";
	}


	@Override
	public int hashCode() {
		return getFile().hashCode();
		}
	
	@Override
	public boolean equals(Object obj) {
		if(obj==null || !(obj instanceof FastQName)) return false;
		if(obj==this ) return true;
		return getFile().equals(((FastQName)obj).getFile());
		}
	
	@Override
	public int compareTo(FastQName o) {
		return getFile().getPath().compareTo(o.getFile().getPath());
		}
	
	public boolean isComplementOf(final FastQName other)
		{
		if(!isValid() || !other.isValid()) return false;
		if(! ((this.isUndetermined() && other.isUndetermined()) || (this.getSample().equals(other.getSample()))) )
			{
			return false;
			}
		
		if(!this.isUndetermined() && !this.getSeqIndex().equals(other.getSeqIndex()))
			{
			return false;
			}
		if(this.getLane()!=other.getLane()) return false;
		if(this.getSplit()!=other.getSplit()) return false;
		if(this.getSeqIndex()!=null && !this.getSeqIndex().equals(other.getSeqIndex())  ) return false;
		if(this.getSide()==other.getSide()) return false;
		return true;
		}
    		
	
	}
