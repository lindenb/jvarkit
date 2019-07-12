package com.github.lindenb.jvarkit.util.bio.structure;

import java.util.List;

import htsjdk.samtools.util.Locatable;

public interface Transcript extends Locatable {
	public Gene getGene();
	public List<Exon> getExons();
    public char getStrand();
    public default int getExonCount() { return this.getExons().size();}
    public default boolean isNegativeStrand() { return getStrand()=='-';}
}
