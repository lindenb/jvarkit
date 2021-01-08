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

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class BcfIteratorBuilder extends VCFIteratorBuilder {

@Override
public VCFIterator open(final String pathOrUrl) throws IOException {
	if(!StringUtils.isBlank(pathOrUrl) && !IOUtil.isUrl(pathOrUrl)) {
		return open(Paths.get(pathOrUrl));
		}
	return super.open(pathOrUrl);
	}
	
@Override
public VCFIterator open(final Path path) throws IOException {
	IOUtil.assertFileIsReadable(path);
	if(BcfToolsUtils.isBcfToolsRequired(path)) {
		final BcfToolsReader br = new BcfToolsReader(path.toString());
		final CloseableIterator<VariantContext> iter = br.iterator();
		final PeekableIterator<VariantContext> peek=new PeekableIterator<>(iter);
		return new VCFIterator() {
			@Override
			public VCFHeader getHeader() {
				return br.getHeader();
				}
			@Override
			public boolean hasNext() {
				return peek.hasNext();
				}
			@Override
			public VariantContext next() {
				return peek.next();
				}
			@Override
			public VariantContext peek() {
				return peek.peek();
				}
			@Override
			public void close() {
				peek.close();
				iter.close();
				try {br.close(); } catch(IOException err) {}
				}
			};
		}
	return super.open(path);
	}


}
