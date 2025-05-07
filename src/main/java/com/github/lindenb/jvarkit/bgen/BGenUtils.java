package com.github.lindenb.jvarkit.bgen;

import java.io.IOException;
import java.nio.charset.Charset;

import htsjdk.samtools.util.BinaryCodec;
import htsjdk.tribble.TribbleException;

public class BGenUtils {
	BGenUtils() {
		
	}
	
	static final Charset ENCODING=Charset.forName("UTF-8");
	static final byte[] BGEN_MAGIC = "bgen".getBytes(ENCODING);
	static final String FILE_SUFFIX = ".bgen";
	
    /** convert a long to int, throws a tribble exception if the number is greater than Integer.MAX_VALUE */
    static int longToUnsignedInt(final long n) {
        if(n>Integer.MAX_VALUE || n<0L) throw new TribbleException("Cannot convert "+n+"L to int");
        return (int)n;
    	}
    /** read the length of a string as unsigned short and read the string */
    static String readStringUInt16(final BinaryCodec binaryCodec) throws IOException {
    	return readString(binaryCodec, binaryCodec.readUShort());
    	}

    	
    /**
     * read a string of length 'len'
     */
   private static String readString(final BinaryCodec binaryCodec, final int len) throws IOException {
        return new String(readNBytes(binaryCodec, len),ENCODING);
    	}
    
    /** fully read 'len' bytes */
   static byte[] readNBytes(final BinaryCodec binaryCodec, final int len) throws IOException {
    	if(len<0) throw new IllegalArgumentException("negative bytes length:" + len);
        final byte[] bytes = new byte[len];
        binaryCodec.readBytes(bytes);
        return bytes;
    	}
}
