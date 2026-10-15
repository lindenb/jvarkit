/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.variant;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFHeader;

/** codec for VariantContext */
public class VariantContextCodec extends AbstractDataCodec<VariantContext> {
	private final VCFHeader header;
    private final VCFCodec vcfDecoder = new VCFCodec();
    private final VCFEncoder vcfEncoder;
    private int capacity = 1000;
	public VariantContextCodec( VCFHeader header) {
		this.header = header;
		this.vcfDecoder.setVCFHeader(header,VCFUtils.getLatestVCFHeaderVersion());
		this.vcfEncoder = new VCFEncoder(this.header, false, true);
		}
	@Override
	public VariantContext decode(final DataInputStream dis) throws IOException {
		final String line = readString(dis);
		this.capacity = Math.max(capacity, line.length());
		return this.vcfDecoder.decode(line);
		}
	
	@Override
	public void encode(final DataOutputStream dos, final VariantContext object) throws IOException {
		final StringBuilder sb=new StringBuilder(this.capacity);
		this.vcfEncoder.write(sb, object);
		writeString(dos, sb.toString());
		}
	
	public VCFHeader getHeader() {
		return header;
		}
	
	@Override
	public VariantContextCodec clone() {
		return new VariantContextCodec(this.header);
		}
	
}
