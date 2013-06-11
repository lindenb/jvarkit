package com.github.lindenb.jvarkit.util.picard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import net.sf.samtools.util.SortingCollection;

public abstract class AbstractDataCodec<T>
	implements SortingCollection.Codec<T>
	{
	private DataInputStream dis=null;
	private DataOutputStream dos=null;
	
		
	public abstract T decode(DataInputStream dis) throws IOException;
	public abstract void encode(DataOutputStream dos,final T object) throws IOException;
	
	@Override
	public abstract AbstractDataCodec<T> clone();
	
	@Override
	public T decode()
		{
		try
			{
			return decode(this.dis);
			}
		catch(IOException err)
			{
			throw new RuntimeException(err);
			}
		}
	
	@Override
	public void encode(final T object)
		{
		try
			{
			encode(this.dos,object);
			}
		catch(IOException err)
			{
			throw new RuntimeException(err);
			}
		}
	
	@Override
	public void setInputStream(InputStream in)
		{
		if(in instanceof DataInputStream)
			{
			this.dis=DataInputStream.class.cast(in);
			}
		else
			{
			this.dis=new DataInputStream(in);
			}
		}
	
	@Override
	public void setOutputStream(OutputStream out)
		{
		if(out instanceof DataOutputStream)
			{
			this.dos=DataOutputStream.class.cast(out);
			}
		else
			{
			this.dos=new DataOutputStream(out);
			}
		}

	
	}
