package com.github.lindenb.jvarkit.pedigree;

import java.util.Set;

/** A Family from a Family */
public interface Family {
/** get family id */
public String getId();
/** get associated pedigree */
public Pedigree getPedigree();
/** return the samples in this family */
public Set<Sample> getSamples();
/** find a sample in this pedigree . returns the sample or null if notr found */
public Sample getSampleById(final String id);
}
