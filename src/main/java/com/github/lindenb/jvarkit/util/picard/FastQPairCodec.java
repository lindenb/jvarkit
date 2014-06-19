package com.github.lindenb.jvarkit.util.picard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.InputStream;
import java.io.OutputStream;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SortingCollection.Codec;

public class FastQPairCodec
implements Codec<FastQPairCodec.PairedFastq>
	{
	private FastQCodec fastqCodec=new FastQCodec();
	private DataInputStream dis;
	private DataOutputStream dos;
	
	public static class PairedFastq
		{
		private FastqRecord fq1;
		private FastqRecord fq2;
		PairedFastq(FastqRecord r1,FastqRecord r2)
			{
			this.fq1=r1;
			this.fq2=r2;
			}
		public FastqRecord getFirst() { return this.fq1;}
		public FastqRecord getSecond() { return this.fq2;}
		}
	
	@Override
	public FastQPairCodec clone()
		{
		return new FastQPairCodec();
		}

	@Override
	public PairedFastq decode()
		{
		FastqRecord fq1=this.fastqCodec.decode();
		if(fq1==null) return null;
		FastqRecord fq2=this.fastqCodec.decode();
		return new PairedFastq(fq1, fq2);	
		}

	@Override
	public void encode(PairedFastq p)
		{
		this.fastqCodec.encode(p.getFirst());
		this.fastqCodec.encode(p.getSecond());
		}

	@Override
	public void setInputStream(InputStream in)
		{
		this.dis=(in instanceof DataInputStream?
				DataInputStream.class.cast(in):
				new DataInputStream(in)
				);
		this.fastqCodec.setInputStream(this.dis);
		}

	@Override
	public void setOutputStream(OutputStream out)
		{
		this.dos=( out instanceof DataOutputStream ?
				DataOutputStream.class.cast(out):
				new DataOutputStream(out)
				);
		this.fastqCodec.setOutputStream(this.dos);
		}
	}
