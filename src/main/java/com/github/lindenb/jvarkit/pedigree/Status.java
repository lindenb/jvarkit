package com.github.lindenb.jvarkit.pedigree;

/** phenotype sample of an individual */
public enum Status {
	missing(-9),unaffected(0),affected(1);

	private final int v;
	Status(int v) { this.v = v;}
	
	public int intValue() { return this.v;}
}
