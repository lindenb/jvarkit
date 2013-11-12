package com.github.lindenb.jvarkit.util.bio;

/**
 * CharSequence of a chromosome 
 * @author lindenb
 *
 */
public interface ChromosomeSequence extends CharSequence
	{
	/** get the chromosome name of that genomic sequence */
	public String getChrom();
	}
