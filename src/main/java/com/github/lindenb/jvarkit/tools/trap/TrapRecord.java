package com.github.lindenb.jvarkit.tools.trap;

import htsjdk.tribble.Feature;

/**
 * A Data Record for the Trap database http://trap-score.org/
 * @author lindenb
 *
 */
public interface TrapRecord extends Feature {
public char getRef();
public char getAlt();
public String getGene();
public float getScore();
}
