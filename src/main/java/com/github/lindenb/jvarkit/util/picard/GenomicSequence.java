package com.github.lindenb.jvarkit.util.picard;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceRecord;

import com.github.lindenb.jvarkit.lang.AbstractCharSequence;

public class GenomicSequence
	extends AbstractCharSequence
	{
	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private SAMSequenceRecord samSequenceRecord;
	
	public GenomicSequence(IndexedFastaSequenceFile indexedFastaSequenceFile ,String chrom)
		{	
		this.indexedFastaSequenceFile=indexedFastaSequenceFile;
		this.samSequenceRecord=this.indexedFastaSequenceFile.getSequenceDictionary().getSequence(chrom);
		}
	
	public SAMSequenceRecord getSAMSequenceRecord()
		{
		return samSequenceRecord;
		}
	
	@Override
	public int hashCode()
		{
		return getSAMSequenceRecord().hashCode();
		}
	
	public String getChrom()
		{
		return getSAMSequenceRecord().getSequenceName();
		}
	
	@Override
	public int length()
		{
		return getSAMSequenceRecord().getSequenceLength();
		}
	
	@Override
	public char charAt(int index0)
		{
		if(index0 >= getSAMSequenceRecord().getSequenceLength())
			{
			throw new IndexOutOfBoundsException("index:"+index0);
			}
		return (char)this.indexedFastaSequenceFile.getSubsequenceAt(getChrom(), index0+1, index0+1).getBases()[0];
		}
	}
