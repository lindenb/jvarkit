package com.github.lindenb.jvarkit.util.jcommander;

import java.util.function.IntSupplier;
import java.util.zip.Deflater;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.StringUtils;

/**
 * Special converter for Zip compression. Bound the values between 0 and 9
 * "best" is interpreted as BEST_COMPRESSION
 * "none" is no compression
 */
public class CompressionConverter  implements IStringConverter<IntSupplier> {
	
	@Override
	public IntSupplier convert(final String s) {
		if(s!=null) {
			if(s.equalsIgnoreCase("best")) return ()->Deflater.BEST_COMPRESSION;
			if(s.equalsIgnoreCase("none")) return ()->Deflater.NO_COMPRESSION;
		}
		if(StringUtils.isBlank(s)) return ()->Deflater.DEFAULT_COMPRESSION;
		final int level = Integer.parseInt(s);
		if(level<Deflater.NO_COMPRESSION) return ()->Deflater.NO_COMPRESSION;
		if(level>Deflater.BEST_COMPRESSION) return ()->Deflater.BEST_COMPRESSION;
			
		return ()->level;
		}
	@Override
	public String toString() {
		return "Compression converter";
		}
	
	
	public static final IntSupplier getDefault() {
		return ()->Deflater.DEFAULT_COMPRESSION;
		}
}