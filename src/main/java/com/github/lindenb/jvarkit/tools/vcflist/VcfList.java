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
package com.github.lindenb.jvarkit.tools.vcflist;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.AbstractList;
import java.util.Collections;
import java.util.List;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Interface used to access a VCF by offset. Static method are used to create a VCF index if needed.
 * @author lindenb
 *
 */
public interface VcfList extends List<VariantContext>,Closeable {
	public VCFHeader getHeader();
	public static VcfList fromFile(final File vcfFile ) throws IOException {
		return fromFile(vcfFile,VcfOffsetsIndexFactory.getDefaultIndexFile(vcfFile));
		}
	public static VcfList fromFile(final File vcfFile,final File indexFile) throws IOException {
		return new DefaultVcfFileList(vcfFile,indexFile);
		}
	public static VcfList from(final VCFHeader header,final List<VariantContext> variants) throws IOException {
		class Tmp extends AbstractList<VariantContext> implements VcfList
			{
			final VCFHeader header;
			final List<VariantContext> variants;
			Tmp(final VCFHeader header,final List<VariantContext> variants)
				{
				this.header=header;
				this.variants = Collections.unmodifiableList(variants);
				}
			@Override
			public VCFHeader getHeader() { return this.header;}
			@Override
			public VariantContext get(int index) {return this.variants.get(index);}
			@Override
			public int size() { return this.variants.size();}
			@Override
			public void close() throws IOException {}
			}
		return new Tmp(header,variants);
		}
	}
