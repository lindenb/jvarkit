package com.github.lindenb.jvarkit.tools.bamstats01;


import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class BamStats01Report implements SAMFileWriter
	{
	private enum Category {
		ALL,
		UNMAPPED,
		MAPPED,
		FAIL_MAPPING_QUALITY,
		DUPLICATE,
		FAIL_VENDOR_QUALITY,
		IN_CAPTURE
		};
	private int qual=0;
	private IntervalTreeMap<? extends Object> bed;
	private Histogram<Category> hist=new Histogram<BamStats01Report.Category>();
	private SAMFileHeader header;
	public BamStats01Report(SAMFileHeader header)
		{
		this.header=header;
		}
	
	public void setBed(IntervalTreeMap<? extends Object> bed)
		{
		this.bed = bed;
		}
	
	@Override
	public void addAlignment(SAMRecord rec)
		{
		if(rec.getReadUnmappedFlag())
			{
			hist.increment(Category.UNMAPPED);
			return;
			}
		hist.increment(Category.MAPPED);
		if(rec.getDuplicateReadFlag())
			{
			hist.increment(Category.DUPLICATE);
			}
		if(rec.getReadFailsVendorQualityCheckFlag())
			{
			hist.increment(Category.FAIL_VENDOR_QUALITY);
			}
		if(qual!=0 && rec.getMappingQuality()<qual)
			{
			hist.increment(Category.FAIL_MAPPING_QUALITY);
			}
		
		if(bed!=null &&
			bed.getOverlapping(new Interval(
					rec.getReferenceName(),
					rec.getAlignmentStart(),
					rec.getAlignmentEnd()
					)).iterator().hasNext()
			)
			{
			hist.increment(Category.IN_CAPTURE);
			}		
		}

	@Override
	public SAMFileHeader getFileHeader()
		{
		return header;
		}

	@Override
	public void close() {
		
		}


}	
