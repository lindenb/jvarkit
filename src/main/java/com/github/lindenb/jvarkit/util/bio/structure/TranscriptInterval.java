package com.github.lindenb.jvarkit.util.bio.structure;

import htsjdk.samtools.util.Locatable;

public interface TranscriptInterval extends Locatable {
public Transcript getTranscript();
public default char getStrand() { return getTranscript().getStrand();}
}
