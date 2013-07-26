package com.github.lindenb.jvarkit.tools.bamstats01;


import net.sf.picard.util.Histogram;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMRecord;

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
