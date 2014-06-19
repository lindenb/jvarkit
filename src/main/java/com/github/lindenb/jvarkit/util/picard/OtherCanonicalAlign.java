package com.github.lindenb.jvarkit.util.picard;

import java.util.List;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

/**  canonical alignments in a chimeric alignment */
public interface OtherCanonicalAlign  extends Comparable<OtherCanonicalAlign>
	{
	@Deprecated
	public String getChrom();
	public String getReferenceName();
	public int getChromIndex();
	@Deprecated
	public int getPos();
	public int getAlignmentStart();
	@Deprecated
	public char getStrand();
	public boolean getReadNegativeStrandFlag();
	public String getCigarString();
	public Cigar getCigar();
	public List<CigarElement> getCigarElements();
	public int getMapQ();
	public int getNM();
	public int getUnclippedStart();
	public int getUnclippedEnd();
	public int getAlignmentEnd();
	}
