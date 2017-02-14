package com.github.lindenb.jvarkit.tools.vcfviewgui;

import java.io.IOException;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Locatable;

public interface NgsFile<HEADERTYPE,ITEMTYPE extends Locatable> {
public HEADERTYPE getHeader();
public SAMSequenceDictionary getSequenceDictionary();
public CloseableIterator<ITEMTYPE> iterator() throws IOException;
public CloseableIterator<ITEMTYPE> iterator(final String contig,final int start,final int end) throws IOException;
public void close();
public String getSource();
public NgsFile<HEADERTYPE,ITEMTYPE> reOpen() throws IOException;
}
