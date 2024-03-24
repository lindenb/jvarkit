package com.github.lindenb.jvarkit.variant.vcf;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.bcf2.BCF2Type;
import htsjdk.variant.bcf2.BCFVersion;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Extract VCF Header from VCF or BCF
 * @author lindenb
 * @see import htsjdk.variant.utils.VCFHeaderReader;
 */

public class VcfHeaderExtractor {
	private VcfHeaderExtractor() {
		
		}
	/** extract VCFHeader from stream. Once decoded, the stream is invalid (uncompressed, bytes skipped, etc..) */
	public static VCFHeader decode(InputStream in0) throws IOException {
		BufferedInputStream in = new BufferedInputStream(IOUtils.mayBeGzippedInputStream(in0),BCFVersion.MAGIC_HEADER_START.length);
        in.mark(2);
		final VCFCodec headerParser = new VCFCodec();

		// read BCF version
        final byte[] magicBytes = new byte[BCFVersion.MAGIC_HEADER_START.length];
        in.read(magicBytes);
        
        if ( !Arrays.equals(magicBytes, BCFVersion.MAGIC_HEADER_START) ) {
        	in.reset();
        	/* looks like a regular VCF */
 	       try( final PositionalBufferedStream bps = new PositionalBufferedStream(in)) {
		        final LineIterator lineIterator = new LineIteratorImpl(new SynchronousLineReader(bps));
		        return (VCFHeader)headerParser.readActualHeader(lineIterator);
		       	}
        	}
        else
	        {
	        // we're a BCF file
	        final int majorByte = in.read();
	        final int minorByte = in.read();
	        new BCFVersion( majorByte, minorByte );
	    	
	        final int headerSizeInBytes = BCF2Type.INT32.read(in);
	
	        if ( headerSizeInBytes <= 0 ) throw new IOException("BCF2 header has invalid length: " + headerSizeInBytes + " must be >= 0");
	        final byte[] headerBytes = new byte[headerSizeInBytes];
	        IOUtils.readFully(in, headerBytes);
	        
	        try( final PositionalBufferedStream bps = new PositionalBufferedStream(new ByteArrayInputStream(headerBytes))) {
		        final LineIterator lineIterator = new LineIteratorImpl(new SynchronousLineReader(bps));
		
		        return (VCFHeader)headerParser.readActualHeader(lineIterator);
		       	}
			}
		}
	
	/** extract VCFHeader from path */
	public static VCFHeader decode(final Path p) throws IOException {
		try(InputStream in = Files.newInputStream(p)) {
			return decode(in);
			}
		}
}
