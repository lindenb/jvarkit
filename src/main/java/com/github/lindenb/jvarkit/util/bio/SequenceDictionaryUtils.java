/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFHeader;

public class SequenceDictionaryUtils {
	
public static boolean isGRCh37(final VCFHeader h) {
return h!=null && isGRCh37(h.getSequenceDictionary());
}

public static boolean isGRCh38(final VCFHeader h) {
return h!=null && isGRCh38(h.getSequenceDictionary());
}

public static boolean isGRCm38(final VCFHeader h) {
return h!=null && isGRCm38(h.getSequenceDictionary());
}
/** return a label for is dictionary or an empty optional */
public static Optional<String> getBuildName(final SAMSequenceDictionary dict) {
	if(dict==null || dict.isEmpty()) return Optional.empty();
	if(isGRCh37(dict)) return Optional.of("GRCh37");
	if(isGRCh38(dict)) return Optional.of("GRCh38");
	if(isGRCm38(dict)) return Optional.of("GRCm38");
	if(isGRCm39(dict)) return Optional.of("GRCm39");
	if(isCanFam3(dict)) return Optional.of("CanFam3");
	if(isCanFam4(dict)) return Optional.of("CanFam4");
	return Optional.empty();
}

private static boolean hasChrom1(final SAMSequenceDictionary dict,int expect) {
	if(dict==null || dict.isEmpty()) return false;
	SAMSequenceRecord rec = dict.getSequence("chr1");
	if(rec==null)  rec = dict.getSequence("1");
	if(rec!=null) {
		if(rec.getSequenceLength()==expect) return true;
		}
	return false;
}

/** test if dict looks like GRCh37  */
public static boolean isGRCh37(final SAMSequenceDictionary dict) {
	return hasChrom1(dict,249_250_621);
	}
/** test if dict looks like GRCh38  */
public static boolean isGRCh38(final SAMSequenceDictionary dict) {
	return hasChrom1(dict,248_956_422);
	}

/** test if dict looks like MusMusculus isGRCm38  */
public static boolean isGRCm38(final SAMSequenceDictionary dict) {
	return hasChrom1(dict,195_471_971);
	}

/** test if dict looks like MusMusculus isGRCm39  */
public static boolean isGRCm39(final SAMSequenceDictionary dict) {
	return hasChrom1(dict,195_154_279);
	}

public static boolean isCanFam3(final VCFHeader h) {
return h!=null && isCanFam3(h.getSequenceDictionary());
}
/** test if dict looks like CanFam3  */
public static boolean isCanFam3(final SAMSequenceDictionary dict) {
	return hasChrom1(dict,122_678_785);
	}

public static boolean isCanFam4(final VCFHeader h) {
return h!=null && isCanFam4(h.getSequenceDictionary());
}
/** test if dict looks like CanFam4  */
public static boolean isCanFam4(final SAMSequenceDictionary dict) {
	return hasChrom1(dict,123_556_469);
	}


/** test if dict looks like a human dict  */
public static boolean isHuman(final VCFHeader h) {
	return isHuman(h.getSequenceDictionary());
	}


/** test if dict looks like a human dict  */
public static boolean isHuman(final SAMSequenceDictionary dict) {
	return isGRCh37(dict) || isGRCh38(dict);
	}


/** test if dict looks like a mouse dict  */
public static boolean isMusMusculus(final VCFHeader h) {
	return isMusMusculus(h.getSequenceDictionary());
	}


/** test if dict looks like a human dict  */
public static boolean isMusMusculus(final SAMSequenceDictionary dict) {
	return isGRCm38(dict) || isGRCm39(dict);
	}


/** extract required from ReferenceSequenceFile */
public static SAMSequenceDictionary extractRequired(final ReferenceSequenceFile faidx) {
	if(faidx==null) throw new IllegalArgumentException("Cannot extract dictionary because ReferenceSequenceFile was not provided.");
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
	if(dict==null || dict.isEmpty())  {
		throw new JvarkitException.VcfDictionaryMissing("<vcf>");
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
			FileExtensions.BAM,
			FileExtensions.SAM,
			FileExtensions.CRAM
			))
			{
			throw new JvarkitException.BamDictionaryMissing(f);
			}
		if(FileExtensions.VCF_LIST.stream().anyMatch(SUFF->f.getFileName().toString().endsWith(SUFF)))
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

/** NO TESTED: TEST THIS. return complement of intervals */
public static List<? extends Locatable> complement(final SAMSequenceDictionary dict,List<? extends Locatable> inputs) {
	final List<Locatable> outputs =new ArrayList<>();
	for(final SAMSequenceRecord ssr:dict.getSequences()) {
		final List<? extends Locatable> L3 = inputs.stream().
				filter(R->R.getContig().equals(ssr.getSequenceName())).
				sorted((A,B)->Integer.compare(A.getStart(), B.getStart())).
				collect(Collectors.toList());
		int pos=1;
		for(int i=0;i< L3.size();i++) {
			final Locatable r=L3.get(i);
			if(!(pos>=r.getStart() && pos<=r.getEnd()))
				{
				outputs.add(new SimpleInterval(ssr.getSequenceName(), pos, r.getStart()-1));
				}
			pos=r.getEnd()+1;
			}
		if(pos< ssr.getSequenceLength()) {
			outputs.add(new SimpleInterval(ssr.getSequenceName(), pos,ssr.getSequenceLength()));
			}
		}
	return outputs;
	}

/** convert a locatable to a queryInterval */
public static Optional<QueryInterval> toQueryInterval(final SAMSequenceDictionary dict,final Locatable loc) {
	if(loc==null) return Optional.empty();
	final SAMSequenceRecord ssr = dict.getSequence(loc.getContig());
	if(ssr==null) return Optional.empty();
	if(loc.getStart()>ssr.getSequenceLength()) return Optional.empty();
	return Optional.of(new QueryInterval(
			ssr.getSequenceIndex(),
			Math.max(1, loc.getStart()),
			Math.min(loc.getEnd(), ssr.getSequenceLength())
			));
	}

}
