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
package com.github.lindenb.jvarkit.samtools;

import java.nio.file.Path;
import java.util.Set;
import java.util.stream.Collectors;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.StringUtil;

/** utilities for samtools */
public class BAMUtils {
/** get first sample name found in header or return default name */
public static String getSampleName(final SAMFileHeader hdr,final String defaultName) {
	return hdr.getReadGroups()
			.stream()
			.map(S->S.getSample())
			.filter(S->!StringUtils.isBlank(S))
			.findAny()
			.orElse(defaultName);
	}
/** get first sample name found in header or return simple name form bamFile */
public static String getSampleName(final SAMFileHeader hdr,final Path bamFile) {
	return getSampleName(hdr,IOUtils.getFilenameWithoutCommonSuffixes(bamFile));
	}
/** return all sample name in a bam file */
public static Set<String> getSampleNames(final SAMFileHeader hdr) {
	return hdr.getReadGroups().stream().map(R->R.getSample())
			.filter(R->!StringUtil.isBlank(R))
			.collect(Collectors.toSet());
	}

/** get mate end position or mate-start, */
public static int getMateEnd(final SAMRecord rec0) {
	if(!rec0.getReadPairedFlag() || rec0.getMateUnmappedFlag()) {
		throw new IllegalStateException("asking mate end but not paired or mate unaligned");
		}
    return rec0.hasAttribute(SAMTag.MC)?
    		SAMUtils.getMateAlignmentEnd(rec0):
    		rec0.getMateAlignmentStart();
    }

}
