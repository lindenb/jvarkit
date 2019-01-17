package com.github.lindenb.jvarkit.pedigree;

import java.util.Collection;

public interface Pedigree {
public Collection<Family> getFamilies();
public default boolean isEmpty() {
	return getFamilies().isEmpty();
	}
public default Family getFamilyById(final String id) {
	return getFamilies().stream().filter(F->F.getId().equals(id)).findAny().orElse(null);
	}
}
