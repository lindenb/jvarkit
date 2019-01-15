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
package com.github.lindenb.jvarkit.util.picard;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;


import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;

public abstract class AbstractDataCodec<T>
	implements SortingCollection.Codec<T>
	{
	private DataInputStream dis=null;
	private DataOutputStream dos=null;
	
		
	public abstract T decode(final  DataInputStream dis) throws IOException;
	public abstract void encode(final DataOutputStream dos,final T object) throws IOException;
	
	@Override
	public abstract AbstractDataCodec<T> clone();
	
	
	
	@Override
	public T decode()
		{
		try
			{
			return decode(this.dis);
			}
		catch(final EOFException err)
			{
			return null;
			}
		catch(final IOException err)
			{
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public void encode(final T object)
		{
		try
			{
			encode(this.dos,object);
			}
		catch(final IOException err)
			{
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public void setInputStream(final  InputStream in)
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
	public void setOutputStream(final OutputStream out)
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
    protected static String readString(final DataInputStream in) throws IOException
    	{
    	return IOUtils.readString(in);
    	}
    
    protected static void writeString(DataOutputStream os,String s) throws IOException
		{
		IOUtils.writeString(os,s);
		}
 

	}
