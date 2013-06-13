package com.github.lindenb.jvarkit.tools.tview;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.SAMRecord;

public abstract class AbstractHandler implements TViewHandler
	{
	public void beginRow() {	 }
	public void endRow() {	 }
	public void beginDocument() {	 }
	public void endDocument() {	 }
	public void beginSAMRecord(final SAMRecord rec) {	 }
	public void endSAMRecord(final SAMRecord rec) {	 }
	public void whiteSpace() {}
	public void deletion() {}
	@Override
	public void base(SAMRecord rec, int readPos, IndexedFastaSequenceFile ref, int refPos) { }
	@Override
	public void beginReferenceSeq() {}
	@Override
	public void endReferenceSeq() {}
	@Override
	public void reference(IndexedFastaSequenceFile ref, String refName,int refPos) {}
	
	protected Character getReferenceBaseAt(IndexedFastaSequenceFile ref, String seqName,int refPos)
		{
		if(ref==null || seqName==null || refPos<0)
			{
			return null;
			}
		ReferenceSequence sub=ref.getSubsequenceAt(seqName,refPos,refPos);
		if(sub==null )return null;
		byte gDNA[]=sub.getBases();
		if(gDNA==null || gDNA.length==0) return null;
		return (char)gDNA[0];
		}

	}
