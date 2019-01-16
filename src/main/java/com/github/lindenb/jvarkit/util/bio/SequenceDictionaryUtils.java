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

import com.github.lindenb.jvarkit.lang.JvarkitException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFHeader;

public class SequenceDictionaryUtils {
	
public static boolean isGRCh37(final VCFHeader h) {
return h!=null && isGRCh37(h.getSequenceDictionary());
}

public static boolean isGRCh38(final VCFHeader h) {
return h!=null && isGRCh38(h.getSequenceDictionary());
}

	
/** test if dict looks like GRCh37  */
public static boolean isGRCh37(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return false;
	SAMSequenceRecord rec = dict.getSequence("chr1");
	if(rec==null)  rec = dict.getSequence("1");
	if(rec!=null) {
		
		if(rec.getSequenceLength()==249250621) return true;
		
		}
	return false;
	}
/** test if dict looks like GRCh38  */
public static boolean isGRCh38(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return false;
	SAMSequenceRecord rec = dict.getSequence("chr1");
	if(rec==null)  rec = dict.getSequence("1");
	if(rec!=null) {
		if(rec.getSequenceLength()==248956422) return true;	
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
public static SAMSequenceDictionary extractRequired(final File f) {
	if(f==null) throw new IllegalArgumentException("Cannot extract dictionary because file was not provided.");
	@SuppressWarnings("deprecation")
	final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(f);
	if(dict==null || dict.isEmpty()) 
		{
		if(f.getName().endsWith(".sam") || f.getName().endsWith(".bam") || f.getName().endsWith(".cram"))
			{
			throw new JvarkitException.BamDictionaryMissing(f);
			}
		if(f.getName().endsWith(".vcf") || f.getName().endsWith(".bcf") || f.getName().endsWith(".vcf.gz"))
			{
			throw new JvarkitException.VcfDictionaryMissing(f);
			}
		throw new JvarkitException.DictionaryMissing(f.getPath());
		}
	return dict;
	}
/** return true if dict has X|chrX AND Y|chrY */
public static boolean hasXY(final SAMSequenceDictionary dict) {
	if(dict==null) return false;
	if(dict.getSequence("X")==null && dict.getSequence("chrX")==null) return false;
	if(dict.getSequence("Y")==null && dict.getSequence("chrY")==null) return false;
	return true;
}
}
