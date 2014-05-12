package com.github.lindenb.jvarkit.tools.tview;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMRecord;

public interface TViewHandler {
	public void beginRow();
	public void endRow();
	public void beginDocument();
	public void endDocument();
	public void beginReferenceSeq();
	public void reference(IndexedFastaSequenceFile ref,String seqName,int refPos);
	public void endReferenceSeq();
	
	public void beginSAMRecord(final SAMRecord rec);
	public void endSAMRecord(final SAMRecord rec);
	public void whiteSpace() ;
	public void deletion();
	public void base(
			SAMRecord rec,int readPos,
			IndexedFastaSequenceFile ref,int refPos
			);

}
