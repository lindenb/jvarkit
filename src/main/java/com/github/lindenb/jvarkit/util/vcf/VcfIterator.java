package com.github.lindenb.jvarkit.util.vcf;

import java.io.Closeable;
import java.io.InputStream;
import java.util.Iterator;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/** net.sf.picard.vcf.VCFIterator deleted from from picard 1.100 */
public class VcfIterator implements Iterator<VariantContext>,Closeable
	{
    private final VCFCodec vcfCodec = new VCFCodec();
    private final VCFHeader vcfHeader;
    private final LineIterator lineIterator;
	public VcfIterator(InputStream vcfStream)
		{
		this( new LineIteratorImpl(LineReaderUtil.fromBufferedStream(vcfStream)));
		}
	
	public VcfIterator(LineIterator r)
		{
		this.lineIterator = r;
	    this.vcfHeader = (VCFHeader) vcfCodec.readActualHeader(lineIterator);
		}
	
	public VCFCodec getCodec()
		{
		return this.vcfCodec;
		}
	
    public VCFHeader getHeader()
    	{
        return this.vcfHeader;
    	}

	
    public VariantContext peek()
    	{
    	return vcfCodec.decode(lineIterator.peek());
		}
    
	@Override
	public boolean hasNext() {
        return lineIterator.hasNext();
		}
	@Override
	public VariantContext next() {
		return vcfCodec.decode(lineIterator.next());
		}
	
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
		}
	
	@Override
	public void close()
		{
		CloserUtil.close(lineIterator);
		}
	}
