package com.github.lindenb.jvarkit.util.picard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import net.sf.picard.fastq.FastqRecord;
import net.sf.samtools.util.SortingCollection.Codec;

public class FastQCodec implements Codec<FastqRecord>
	{
	private DataInputStream dis;
	private DataOutputStream dos;
	private byte buffer[]=new byte[100];
	
	@Override
	public FastQCodec clone()
		{
		return new FastQCodec();
		}

	private void writeAscii(final String s)
		throws IOException
		{
		byte a[]=s.getBytes();
		this.dos.writeInt(a.length);
		this.dos.write(a);
		}
	
	private String readAscii()
		throws IOException
		{
		int n=this.dis.readInt();
		if(this.buffer.length<n) this.buffer=new byte[n*2];
		this.dis.read(this.buffer,0,n);
		return new String(this.buffer,0,n);
		}

	
	
	@Override
	public FastqRecord decode()
		{
		try
			{
			if(this.dis.available()==0) return null;
			String L1=readAscii();
			String L2=readAscii();
			String L3=readAscii();
			String L4=readAscii();
			return new FastqRecord(L1,L2,L3,L4);
			}
		catch (IOException e)
			{
			throw new RuntimeException(e);
			}
		}

	@Override
	public void encode(FastqRecord p)
		{
		try
			{
			writeAscii(p.getReadHeader());
			writeAscii(p.getReadString());
			writeAscii(p.getBaseQualityHeader());
			writeAscii(p.getBaseQualityString());
			}
		catch (IOException e)
			{
			throw new RuntimeException(e);
			}
		}

	@Override
	public void setInputStream(InputStream in)
		{
		this.dis=(in instanceof DataInputStream?
				DataInputStream.class.cast(in):
				new DataInputStream(in)
				);
		}

	@Override
	public void setOutputStream(OutputStream out)
		{
		this.dos=( out instanceof DataOutputStream ?
				DataOutputStream.class.cast(out):
				new DataOutputStream(out)
				);
		}
	}
