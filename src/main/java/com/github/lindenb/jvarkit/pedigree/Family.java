package com.github.lindenb.jvarkit.pedigree;

import java.util.Set;

public interface Family {
public String getId();
public Set<Sample> getSamples();
public default Sample getSampleById(final String id) {
	return getSamples().stream().filter(S->S.getId().equals(id)).findAny().orElse(null);
	}
}
