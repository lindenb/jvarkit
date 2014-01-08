package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.OutputStream;
/** output stream that doesn't print anything */
public class NullOuputStream extends OutputStream
	{
	private boolean _closed=false;
	private long count=0L;
	public NullOuputStream()
		{
		}
	@Override
	public void write(int c)// throws IOException
		{
		if(!isClosed() && c>=0) count++;
		}
	@Override
	public void write(byte[] b)// throws IOException
		{
		if(!isClosed()) count+=b.length;
		}
	@Override
	public void write(byte[] b, int off, int len)// throws IOException
		{
		if(!isClosed()) count+=len;
		}
	@Override
	public void close() throws IOException
		{
		_closed=true;
		}
	public long getByteWrittenCount()
		{
		return count;
		}
	public boolean isClosed()
		{
		return _closed;
		}
	@Override
	public String toString()
		{
		return "NullOuputStream: closed:"+ isClosed()+" N="+getByteWrittenCount();
		}
	}
