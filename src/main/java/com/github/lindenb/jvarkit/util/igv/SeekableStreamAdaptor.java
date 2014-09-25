package com.github.lindenb.jvarkit.util.igv;
import java.io.IOException;

import htsjdk.samtools.seekablestream.SeekableStream;
/**
 * Motivation: use SeekableStream from htsjdk.samtools.seekablestream
 * to org.broad.tribble.util.SeekableStream 
 */
public class SeekableStreamAdaptor
extends org.broad.tribble.util.SeekableStream
{
private SeekableStream delegate;

public SeekableStreamAdaptor(SeekableStream delegate)
	{
	this.delegate=delegate;
	}

private SeekableStream getDelegate()
	{
	return this.delegate;
	}

public int available() throws IOException {
	return getDelegate().available();
}

public void close() throws IOException {
	getDelegate().close();
}

public boolean eof() throws IOException {
	return getDelegate().eof();
}


@Override
public int hashCode() {
	return getDelegate().hashCode();
}

@Override
public long length() {
	return getDelegate().length();
}

@Override
public void mark(int readlimit) {
	getDelegate().mark(readlimit);
}

@Override
public boolean markSupported() {
	return getDelegate().markSupported();
}

@Override
public long position() throws IOException {
	return getDelegate().position();
}

@Override
public int read() throws IOException {
	return getDelegate().read();
}

@Override
public int read(byte[] arg0, int arg1, int arg2) throws IOException {
	return getDelegate().read(arg0, arg1, arg2);
}

@Override
public int read(byte[] b) throws IOException {
	return getDelegate().read(b);
}

@Override
public void readFully(byte[] arg0) throws IOException {
	getDelegate().readFully(arg0);
}

@Override
public void reset() throws IOException {
	getDelegate().reset();
}

@Override
public void seek(long arg0) throws IOException {
	getDelegate().seek(arg0);
}

@Override
public long skip(long n) throws IOException {
	return getDelegate().skip(n);
}

}
