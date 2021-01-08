/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.InputStream;
import java.io.PushbackInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.GZIPInputStream;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.IOUtil;

/** uncompress delegate IF it's look like a Gzipped file */
public class MayBeGzipInputStream extends InputStream {

private final InputStream decompressor;	

public MayBeGzipInputStream(final Path path) throws IOException {
	this(path,512);
	}
public MayBeGzipInputStream(final Path path,int buffSize) throws IOException {
	IOUtil.assertFileIsReadable(path);
	final byte [] signature = new byte[2];
	final int n;
	try(InputStream fin = Files.newInputStream(path)) {
		n= fin.read( signature ); //read the signature
		}
	if(n!=2 || !IOUtils.isGZipCompressed(signature)) {
    	this.decompressor = Files.newInputStream(path) ;
     	}
     else
     	{
    	this.decompressor =  new GZIPInputStream( Files.newInputStream(path),buffSize );
     	}
	}

public MayBeGzipInputStream(final InputStream delegate) throws IOException {
	this(delegate,512);
	}

public MayBeGzipInputStream(final InputStream delegate,int buffSize) throws IOException {
	 if(delegate==null) throw new IllegalArgumentException("input is null");
	 if(delegate instanceof GZIPInputStream) {
		decompressor = delegate;
	 	}
	 else if(delegate instanceof BlockCompressedInputStream) {
		decompressor = delegate;
	 	}
	 else
		 {	
		 final byte [] signature = new byte[2];
	     final PushbackInputStream pb = new PushbackInputStream(delegate, signature.length ); //we need a pushbackstream to look ahead
	     final int n= pb.read( signature ); //read the signature
	     pb.unread( signature,0,n);
	     if(n!=2 || !IOUtils.isGZipCompressed(signature)) {
	    	this.decompressor = pb;
	     	}
	     else
	     	{
	    	this.decompressor =  new GZIPInputStream( pb,buffSize );
	     	}
		}
	}


@Override
public int read() throws IOException {
	return this.decompressor.read();
	}

@Override
	public int read(byte[] b, int off, int len) throws IOException {
	return this.decompressor.read(b, off, len);
	}

 @Override
	public void close() throws IOException {
	this.decompressor.close();
	}

}
