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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;

/** Factory creating VCFReader */
public abstract class VCFReaderFactory {
	private boolean requireIndex = true; /* because default is true in VCFFileReader */
	
	/** initialize with requireIndex = true */
	protected VCFReaderFactory() {	
		}
	
	public static VCFReaderFactory makeDefault() {
		return new VCFReaderFactory() {
			};
		}
	
	public VCFReaderFactory setRequireIndex(boolean b) {
		this.requireIndex = b;
		return this;
		}
	
	public boolean isRequireIndex() {
		return requireIndex;
		}
		
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public VCFReader open(final String pathOrUrl) {
		return open(pathOrUrl,isRequireIndex());
	}
	
	/** open new VCFReader */
	public VCFReader open(final String pathOrUrl,boolean requireIndex) {
		if(pathOrUrl.startsWith("file://")) {
			return open(pathOrUrl.substring(7),requireIndex);
			}
		if(IOUtil.isUrl(pathOrUrl) ) {
			if(requireIndex) {
				return new TabixVcfReader(pathOrUrl);
				}
			else
				{	
				try {
					return new SimpleVcfIteratorWrapper(  new BcfIteratorBuilder().open(pathOrUrl));
					}
				catch(IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			}

		return open(Paths.get(pathOrUrl), requireIndex);
		}

	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public VCFReader open(final Path p) {
		return open(p,isRequireIndex());
		}
	
	/** open new VCFReader */
	public VCFReader open(final Path path,boolean requireIndex) {
		if(BcfToolsUtils.isBcfToolsRequired(path)) {
			return new BcfToolsReader(path.toString());
			}
		
		return new VCFFileReader(path, requireIndex);
		}
	
	/** open new VCFReader */
	public final VCFReader open(final File path,boolean requireIndex) {
		return open(path.toPath(),requireIndex);
		}
	
	/** open new VCFReader with default {@link #isRequireIndex()} */
	public final VCFReader open(final File p) {
		return open(p,isRequireIndex());
		}

	private static class SimpleVcfIteratorWrapper implements VCFReader {
		private final VCFIterator iter;
		SimpleVcfIteratorWrapper(VCFIterator iter) {
			this.iter = iter;
		 	}
		@Override
		public VCFHeader getHeader() {
			return this.iter.getHeader();
			}
		@Override
		public boolean isQueryable() {
			return false;
			}
		@Override
		public CloseableIterator<VariantContext> iterator() {
			return this.iter;
			}
		@Override
		public CloseableIterator<VariantContext> query(String chrom, int start, int end) {
			throw new IllegalStateException("Cannot query a VCF iterator because index was not required");
			}
		@Override
		public void close() throws IOException {
			this.iter.close();
			}
		}
	
}
