package com.github.lindenb.jvarkit.util.picard;

import java.util.List;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;

/**  canonical alignments in a chimeric alignment */
public interface OtherCanonicalAlign  extends Comparable<OtherCanonicalAlign>
	{
	public String getChrom();
	public int getChromIndex();
	public int getPos();
	public char getStrand();
	public String getCigarString();
	public Cigar getCigar();
	public List<CigarElement> getCigarElements();
	public int getMapQ();
	public int getNM();
	}
