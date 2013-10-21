package com.github.lindenb.jvarkit.util.picard;

import java.util.List;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;


public interface XPAlign  extends Comparable<XPAlign>
	{
	public String getChrom();
	public int getChromIndex();
	public int getPos();
	public char getStrand();
	public String getCigarString();
	public Cigar getCigar();
	public List<CigarElement> getCigarElements();
	}
