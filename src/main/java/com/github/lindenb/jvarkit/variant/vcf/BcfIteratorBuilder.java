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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.variant.bcf.BCFIterator;

import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.bcf2.BCFVersion;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

/** custom Builder making VCFIterators */
public class BcfIteratorBuilder extends VCFIteratorBuilder {

	@Override
	public VCFIterator open(final InputStream in) throws IOException {
        if (in == null) {
            throw new IllegalArgumentException("input stream is null");
            }
        // wrap the input stream into a BufferedInputStream to reset/read a BCFHeader or a GZIP
        // buffer must be large enough to contain the BCF header and/or GZIP signature
        BufferedInputStream  bufferedinput = new BufferedInputStream(in, Math.max(BCF2Codec.SIZEOF_BCF_HEADER, IOUtil.GZIP_HEADER_READ_LENGTH));
        // test for gzipped inputstream
        if(IOUtil.isGZIPInputStream(bufferedinput)) {
            // this is a gzipped input stream, wrap it into GZIPInputStream
            // and re-wrap it into BufferedInputStream so we can test for the BCF header
            bufferedinput = new BufferedInputStream(new GZIPInputStream(bufferedinput), BCF2Codec.SIZEOF_BCF_HEADER);
        }

        // try to read a BCF header
        final BCFVersion bcfVersion = BCF2Codec.tryReadBCFVersion(bufferedinput);

        if (bcfVersion != null) {
            //this is BCF
            return BCFIterator.open(bufferedinput);
        } else {
            //this is VCF
            return new VCFReaderIterator2(bufferedinput);
        }
    }	
	
	
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
	if(path.getFileName().getFileName().toString().endsWith(FileExtensions.BCF)) {
		return BCFIterator.open(path);
		}
	return super.open(path);
	}



/** implementation of VCFIterator, reading VCF, sadly super.VCFReaderIterator is private
 * TODO : remove this class when it  is public is htsjdk
 */
private static class VCFReaderIterator2
        extends AbstractIterator<VariantContext>
        implements VCFIterator {
    /** delegate input stream */
    private final InputStream inputStream;
    /** VCF codec */
    private final VCFCodec codec = new VCFCodec();
    /** VCF header */
    private final VCFHeader vcfHeader;
    /** Iterator over the lines of the VCF */
    private final LineIterator lineIterator;

    VCFReaderIterator2(final InputStream inputStream) {
        this.inputStream = inputStream;
        this.lineIterator = this.codec.makeSourceFromStream(this.inputStream);
        this.vcfHeader = (VCFHeader) this.codec.readActualHeader(this.lineIterator);
    }

    @Override
    public VCFHeader getHeader() {
        return this.vcfHeader;
    }

    @Override
    protected VariantContext advance() {
        return this.lineIterator.hasNext() ? this.codec.decode(this.lineIterator.next()) : null;
    }

    @Override
    public void close() {
    	try{this.inputStream.close();} catch(Throwable err) {}
    	}
	}
}
