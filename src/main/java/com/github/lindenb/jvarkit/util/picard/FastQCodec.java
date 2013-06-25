package com.github.lindenb.jvarkit.util.picard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import net.sf.picard.fastq.FastqRecord;

public class FastQCodec extends AbstractDataCodec<FastqRecord>
	{
	private byte buffer[]=new byte[100];
	
	@Override
	public FastQCodec clone()
		{
		return new FastQCodec();
		}

	private void writeAscii(DataOutputStream dos,final String s)
		throws IOException
		{
		byte a[]=s.getBytes();
		dos.writeInt(a.length);
		dos.write(a);
		}
	
	private String readAscii(DataInputStream dis)
		throws IOException
		{
		int n=dis.readInt();
		if(this.buffer.length<n) this.buffer=new byte[n*2];
		dis.read(this.buffer,0,n);
		return new String(this.buffer,0,n);
		}

	@Override
	public FastqRecord decode(DataInputStream dis) throws IOException
		{
		String L1;
		try
			{
			 L1=readAscii(dis);
			}
		catch(IOException err)
			{
				return null;
			}
		String L2=readAscii(dis);
		String L3=readAscii(dis);
		String L4=readAscii(dis);
		return new FastqRecord(L1,L2,L3,L4);
		}
	
	@Override
	public void encode(DataOutputStream dos, FastqRecord p)
			throws IOException
		{
		writeAscii(dos,p.getReadHeader());
		writeAscii(dos,p.getReadString());
		writeAscii(dos,p.getBaseQualityHeader());
		writeAscii(dos,p.getBaseQualityString());
		}


	}
