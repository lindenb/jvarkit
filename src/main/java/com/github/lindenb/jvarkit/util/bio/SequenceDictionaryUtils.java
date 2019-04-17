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
package com.github.lindenb.jvarkit.util.bio;

import java.io.File;
import java.nio.file.Path;
import java.util.Optional;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFHeader;

public class SequenceDictionaryUtils {
	
public static boolean isGRCh37(final VCFHeader h) {
return h!=null && isGRCh37(h.getSequenceDictionary());
}

public static boolean isGRCh38(final VCFHeader h) {
return h!=null && isGRCh38(h.getSequenceDictionary());
}

/** return a label for is dictionary or an empty optional */
public static Optional<String> getBuildName(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return Optional.empty();
	if(isGRCh37(dict)) return Optional.of("GRCh37");
	if(isGRCh38(dict)) return Optional.of("GRCh38");
	return Optional.empty();
}

/** test if dict looks like GRCh37  */
public static boolean isGRCh37(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return false;
	SAMSequenceRecord rec = dict.getSequence("chr1");
	if(rec==null)  rec = dict.getSequence("1");
	if(rec!=null) {
		
		if(rec.getSequenceLength()==249_250_621) return true;
		
		}
	return false;
	}
/** test if dict looks like GRCh38  */
public static boolean isGRCh38(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return false;
	SAMSequenceRecord rec = dict.getSequence("chr1");
	if(rec==null)  rec = dict.getSequence("1");
	if(rec!=null) {
		if(rec.getSequenceLength()==248_956_422) return true;	
		}
	return false;
	}

/** test if dict looks like a human dict  */
public static boolean isHuman(final VCFHeader h) {
	return isHuman(h.getSequenceDictionary());
	}


/** test if dict looks like a human dict  */
public static boolean isHuman(final SAMSequenceDictionary dict) {
	return isGRCh37(dict) || isGRCh38(dict);
	}

/** extract required from IndexedFastaSequenceFile */
public static SAMSequenceDictionary extractRequired(final IndexedFastaSequenceFile faidx) {
	if(faidx==null) throw new IllegalArgumentException("Cannot extract dictionary because IndexedFastaSequenceFile was not provided.");
	final SAMSequenceDictionary dict = faidx.getSequenceDictionary();
	if(dict==null || dict.isEmpty()) 
		{
		throw new JvarkitException.FastaDictionaryMissing(faidx.toString());
		}
	return dict;
	}


/** extract required SAMSequenceDictionary */
public static SAMSequenceDictionary extractRequired(final SAMFileHeader h) {
	if(h==null) throw new IllegalArgumentException("Cannot extract dictionary because SAM header was not provided.");
	final SAMSequenceDictionary dict = h.getSequenceDictionary();
	if(dict==null || dict.isEmpty()) 
		{
		throw new JvarkitException.BamDictionaryMissing("<bam>");
		}
	return dict;
	}

/** extract required SAMSequenceDictionary */
public static SAMSequenceDictionary extractRequired(final VCFHeader h) {
	if(h==null) throw new IllegalArgumentException("Cannot extract dictionary because VCF header was not provided.");
	final SAMSequenceDictionary dict = h.getSequenceDictionary();
	if(dict==null || dict.isEmpty()) 
		{
		throw new JvarkitException.BamDictionaryMissing("<vcf>");
		}
	return dict;
	}

/** extract required SAMSequenceDictionary */
public static SAMSequenceDictionary extractRequired(final Path f) {
	if(f==null) throw new IllegalArgumentException("Cannot extract dictionary because file was not provided.");
	final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(f);
	if(dict==null || dict.isEmpty()) 
		{
		if(StringUtils.endsWith(f.getFileName().toString(),
			BamFileIoUtils.BAM_FILE_EXTENSION,
			IOUtil.SAM_FILE_EXTENSION,
			CramIO.CRAM_FILE_EXTENSION
			))
			{
			throw new JvarkitException.BamDictionaryMissing(f);
			}
		if(StringUtils.endsWith(f.getFileName().toString(),
				IOUtil.VCF_FILE_EXTENSION,
				IOUtil.COMPRESSED_VCF_FILE_EXTENSION,
				IOUtil.BCF_FILE_EXTENSION)
				)
			{
			throw new JvarkitException.VcfDictionaryMissing(f.toString());
			}
		throw new JvarkitException.DictionaryMissing(f.toString());
		}
	return dict;
	}
/** return true if dict has X|chrX AND Y|chrY */


/** extract required SAMSequenceDictionary */
public static SAMSequenceDictionary extractRequired(final File f) {
	if(f==null) throw new IllegalArgumentException("Cannot extract dictionary because file was not provided.");
	return extractRequired(f.toPath());
	}

/** return true if dict has X|chrX AND Y|chrY */
public static boolean hasXY(final SAMSequenceDictionary dict) {
	if(dict==null) return false;
	if(dict.getSequence("X")==null && dict.getSequence("chrX")==null) return false;
	if(dict.getSequence("Y")==null && dict.getSequence("chrY")==null) return false;
	return true;
}
}
