package com.github.lindenb.jvarkit.util.illumina;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;

import net.sf.picard.fastq.FastqReader;

/** FastQName */
public class FastQName
	implements Comparable<FastQName>
	{
	private static final String SUFFIX=".fastq.gz";
	private File file;
	private boolean is_valid=false;
	private String seqIndex;
	private int lane;
	private int side=-1;
	private String sample=null;
	private int split=-1;
	private boolean undetermined=false;
	private Long nReads=null;
	
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
	
	public int getSide()
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
	
	//SAMPLENAME_GATCAG_L007_R2_001.fastq.gz
	//lane7_Undetermined_L007_R2_008.fastq.gz
	int i0=f.getName().lastIndexOf('_');
	if(i0==-1)
		{
		return fq;
		}
	
	String splitS=f.getName().substring(i0+1);
	fq.split=Integer.parseInt(splitS.substring(0,splitS.length()-SUFFIX.length()));
	
	int i1=f.getName().lastIndexOf('_',i0-1);
	if(i1==-1)
		{
		return fq;
		}
	
	String sideStr=f.getName().substring(i1+1,i0);
	if(sideStr.equals("R1"))
		{
		fq.side=0;
		}
	else if(sideStr.equals("R2"))
		{
		fq.side=1;
		}
	else
		{
		return fq;
		}
	
	int i2=f.getName().lastIndexOf('_',i1-1);
	if(i2==-1)
		{
		return fq;
		}
	String laneStr=f.getName().substring(i2+1,i1);
	if(!laneStr.startsWith("L"))
		{
		return fq;
		}
	fq.lane=Integer.parseInt(laneStr.substring(1));
	
	if(f.getName().substring(0,i2).toLowerCase().endsWith("undetermined"))
		{
		fq.is_valid=true;
		fq.undetermined=true;
		return fq;
		}
	
	int i3=f.getName().lastIndexOf('_',i2-1);
	if(i3==-1)
		{
		return fq;
		}
	fq.seqIndex=f.getName().substring(i3+1,i2);
	
	fq.sample=f.getName().substring(0, i3);
	
	
	
	
		
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
				r=new FastqReader(this.getFile());
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
		if(! this.getSeqIndex().equals(other.getSeqIndex())  ) return false;
		if(this.getSide()==other.getSide()) return false;
		return true;
		}
    		
	
	}
