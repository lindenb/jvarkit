package com.github.lindenb.jvarkit.util.bio.structure;

import java.util.List;
import java.util.Map;

import htsjdk.samtools.util.Locatable;

public interface Gene extends Locatable {
public List<Transcript> getTranscripts();
public Map<String, String> getProperties();
}
