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

import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMUtils;

public class SAMRecordRead implements Read {
private final SAMRecord delegate;
public SAMRecordRead(final SAMRecord delegate) {
	this.delegate = Objects.requireNonNull(delegate);
	}
public SAMRecord getDelegate() {
	return delegate;
	}
@Override
public String getContig() {
	return getDelegate().getContig();
	}
@Override
public int getStart() {
	return getDelegate().getStart();
	}
@Override
public int getEnd() {
	return getDelegate().getEnd();
	}
@Override
public String getName() {
	return getDelegate().getReadName();
	}
@Override
public boolean isDuplicate() {
	return getDelegate().getDuplicateReadFlag();
	}
@Override
public int getMappingQuality() {
	return getDelegate().getMappingQuality();
	}
@Override
public boolean isSecondary() {
	return getDelegate().isSecondaryAlignment();
	}
@Override
public boolean isUnmapped() {
	return getDelegate().getReadUnmappedFlag();
	}
@Override
public boolean isFailingVendorQualityCheck() {
	return getDelegate().getReadFailsVendorQualityCheckFlag();
	}
@Override
public boolean isSupplementary() {
	return getDelegate().getSupplementaryAlignmentFlag();
	}
@Override
public boolean isPairedRead() {
	return getDelegate().getReadPairedFlag();
	}
@Override
public Cigar getCigar() {
	return getDelegate().getCigar();
}
@Override
public boolean isMateUnmapped() {
	return getDelegate().getMateUnmappedFlag();
}
@Override
public int getReferenceIndex() {
	return getDelegate().getReferenceIndex().intValue();
}
@Override
public int getMateReferenceIndex() {
	return getDelegate().getMateReferenceIndex();
}
@Override
public String getMateReferenceName() {
	return getDelegate().getMateReferenceName();
}
@Override
public int getMateAlignmentStart() {
	return getDelegate().getMateAlignmentStart();
}
@Override
public int getMateAlignmentEnd() {
	return BAMUtils.getMateEnd(getDelegate());
}
@Override
public List<Read> getOtherCanonicalAlignments() {
	return SAMUtils.getOtherCanonicalAlignments(getDelegate()).stream().map(R->new SAMRecordRead(R)).collect(Collectors.toList());
	}

@Override
public boolean isPlusStrand() {
	return !isMinusStrand();
}
@Override
public boolean isMinusStrand() {
	return getDelegate().getReadNegativeStrandFlag();
}
@Override
public boolean isMatePlusStrand() {
	return !isMateMinusStrand();
}
@Override
public boolean isMateMinusStrand() {
	return getDelegate().getMateNegativeStrandFlag();
}
}
