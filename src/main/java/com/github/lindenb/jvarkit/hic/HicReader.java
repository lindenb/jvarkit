package com.github.lindenb.jvarkit.hic;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;

public interface HicReader extends AutoCloseable {
	
/** get source of this reader (path, url...) or null */
public Object getSource();	

/** get dictionary */
public SAMSequenceDictionary getDictionary();

/** get genome build */
public String getBuild();

/** get attributes */
public Map<String,String> getAttributes();

/** get version of hic format */
public int getVersion();

/** get the base pair resolutions */
public Set<Integer> getBasePairResolutions();

/** get the fragment resolutions */
public Set<Integer> getFragmentResolutions();

public Iterator<QueryResult> query(
		final Locatable interval1,
		final Locatable interval2,
		final Normalization norm,
		final int binsize, 
		final Unit unit
		);

public static HicReader open(final Path path) throws IOException {
	return new HicReaderImpl(path);
	}
}
