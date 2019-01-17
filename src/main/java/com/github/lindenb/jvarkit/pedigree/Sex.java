package com.github.lindenb.jvarkit.pedigree;

/** sex of an individual */
public enum Sex {
	male(1),female(2),unknown(0);
	private final int v;
	Sex(int v) {
		this.v = v;
		}
	
	public int intValue() { return this.v;}
	}
