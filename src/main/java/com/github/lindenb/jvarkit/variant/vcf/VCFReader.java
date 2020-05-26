/*
The MIT License (MIT)

Copyright (c) 2020 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.File;
import java.io.*;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.vcf.*;

/** DELETE WHEN HTSJDK is upgraded */
public interface VCFReader extends Closeable {
	

	public VCFHeader getHeader();
	public CloseableIterator<VariantContext> query(String contig, int start, int end);
	public CloseableIterator<VariantContext> iterator();

	public static VCFReader open(final Path path) {
		try {
			if(path.toString().endsWith(".bcf")) {
				return new BcfToolsReader(path.toString());
			} else {
				return new VCFFileReaderWrapper(path);
			}
		  }
		catch(Exception err) {
			throw new RuntimeException(err);
			}
		}

	static class VCFFileReaderWrapper implements VCFReader {
		private final VCFFileReader delegate;
		VCFFileReaderWrapper(final Path path) throws IOException {
			this.delegate = new VCFFileReader(path,true);
			}
		@Override
		public VCFHeader getHeader() { return this.delegate.getFileHeader();}
		@Override
		public CloseableIterator<VariantContext> iterator() { return this.delegate.iterator();}
		@Override
		public CloseableIterator<VariantContext> query(String contig, int start, int end) { return this.delegate.query(contig,start,end);}
		@Override
		public void close() throws IOException { try{ this.delegate.close();} catch(Exception err){}}
	}

	
}

