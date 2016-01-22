package com.github.lindenb.jvarkit.util.vcf;

import java.io.InputStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReaderUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.AbstractVCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/** net.sf.picard.vcf.VCFIterator deleted from from picard 1.100 */
public class VcfIteratorImpl implements VcfIterator
	{
	/** associated VCF codec */
    private final AbstractVCFCodec vcfCodec;
	/** associated VCF header */
    private final VCFHeader vcfHeader;
	/** associated line iterator */
    private final LineIterator lineIterator;
    
	public VcfIteratorImpl(InputStream vcfStream)
		{
		this( new LineIteratorImpl(LineReaderUtil.fromBufferedStream(vcfStream)));
		}
	
	public VcfIteratorImpl(LineIterator r)
		{
		this.lineIterator = r;
	   // this.vcfHeader = (VCFHeader) vcfCodec.readActualHeader(lineIterator);
	    
	    VCFUtils.CodecAndHeader cah=VCFUtils.parseHeader(r);
	    this.vcfHeader=cah.header;
	    this.vcfCodec=cah.codec;
	    
		}
	
	@Override
	public AbstractVCFCodec getCodec()
		{
		return this.vcfCodec;
		}
	
	@Override
    public VCFHeader getHeader()
    	{
        return this.vcfHeader;
    	}

	@Override
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
	@Override
	public String toString() {
		return "VCF Iterator. Codec: "+vcfCodec;
		}
	}
