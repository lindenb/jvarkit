package com.github.lindenb.jvarkit.util.vcf;

import java.io.Closeable;
import java.io.InputStream;
import java.util.Iterator;
import net.sf.samtools.util.CloserUtil;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReaderUtil;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.vcf.VCFCodec;
import org.broadinstitute.variant.vcf.VCFHeader;

/** net.sf.picard.vcf.VCFIterator deleted from from picard 1.100 */
public class VcfIterator implements Iterator<VariantContext>,Closeable
	{
    private final VCFCodec vcfCodec = new VCFCodec();
    private final VCFHeader vcfHeader;
    private final LineIterator lineIterator;
	public VcfIterator(InputStream vcfStream)
		{
		this.lineIterator = new LineIteratorImpl(LineReaderUtil.fromBufferedStream(vcfStream));
	    this.vcfHeader = (VCFHeader) vcfCodec.readActualHeader(lineIterator);
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
