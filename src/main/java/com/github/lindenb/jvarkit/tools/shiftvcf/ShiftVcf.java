package com.github.lindenb.jvarkit.tools.shiftvcf;

import java.nio.file.Path;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC



shift coordinates of a VCF.


## Usage 1

shift a vcf using a using ROI.bed

```
bcftools view --regions-file ROI.bed source.vcf.gz |\
	java -jar dist/jvarkit.jar shiftvcf --bed ROI.bed
```


## Usage 2 

Shift the VCF back to original coordinates


```
java -jar dist/jvarkit.jar shiftvcf --R2 source.vcf.gz shifted.vcf

````




END_DOC
*/
@Program(name="shiftvcf",
description="shit all coordinates of a VCF",
keywords={"vcf","bed"},
creationDate="20241002",
modificationDate="20241002",
jvarkit_amalgamion = true
)
public class ShiftVcf extends OnePassVcfLauncher {
	private static final Logger LOG=Logger.build(ShiftVcf.class).make(); 
	@Parameter(names={"--bed"},description="Bed Path. Extract Variants overlaping this BED. Or use -R2.")
	private Path bedPath = null;
	@Parameter(names={"-R2","--destination-reference"},description="Original fasta reference. We shift the VCF back to this reference. Required without --bed." + DICTIONARY_SOURCE)
	private Path destinationRefPath = null;
	
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	private static Locatable contig2interval(final String s) {
		final int colon=s.lastIndexOf(':');
		if(colon==-1) throw new IllegalArgumentException("cannot find ':' in "+s);
		int hyphen=s.lastIndexOf('-');
		if(hyphen==-1 || hyphen < colon) throw new IllegalArgumentException("cannot find '-' in "+s);
		final String contig = s.substring(0,colon);
		final int start = Integer.parseInt(s.substring(colon+1,hyphen));
		final int end = Integer.parseInt(s.substring(hyphen+1));
		if(start>end) throw new IllegalArgumentException("start> end in "+s);
		return new SimpleInterval(contig, start, end);
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		if(this.bedPath==null) {
			return doWorkWithoutBed(inputName,iterin,out);
			}
		else
			{
			return doWorkWithBed(inputName,iterin,out);
			}
		}
	
	
	
	private int doWorkWithoutBed(String inputName, VCFIterator iterin, VariantContextWriter out) {
		try {
			if(this.destinationRefPath==null) {
				LOG.error("destination reference path is required.");
				return -1;
				}
			final VCFHeader headerIn  = iterin.getHeader();

			final SAMSequenceDictionary destDict = new SequenceDictionaryExtractor().
					extractRequiredDictionary(this.destinationRefPath);
			
			
			final VCFHeader headerOut = new VCFHeader(headerIn);
			headerOut.setSequenceDictionary(destDict);
			JVarkitVersion.getInstance().addMetaData(this, headerOut);
			out.writeHeader(headerOut);
			while(iterin.hasNext()) {
				final VariantContext ctx0 = iterin.next();
				
				final Locatable loc= contig2interval(ctx0.getContig());
				if(destDict.getSequence(loc.getContig())==null) throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(),destDict);

				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx0);
				vcb.chr(loc.getContig());
				vcb.start(loc.getStart() + ctx0.getStart() -1);
				final int chromEnd= loc.getStart() + ctx0.getEnd() -1;
				vcb.stop(chromEnd);
				if(ctx0.hasAttribute(VCFConstants.END_KEY)) {
					vcb.attribute(VCFConstants.END_KEY,chromEnd);
					}
				out.add(vcb.make());
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	
	private int doWorkWithBed(String inputName, VCFIterator iterin, VariantContextWriter out) {
		try {		
			final VCFHeader headerIn  = iterin.getHeader();
			final SAMSequenceDictionary dict0 =  SequenceDictionaryUtils.extractRequired(headerIn);
			
			final IntervalTreeMap<Interval> intervalTreeMap=new IntervalTreeMap<>();
			try(BedLineReader rd = new BedLineReader(this.bedPath)) {
				while(rd.hasNext()) {
					final BedLine bedrec = rd.next();
					final Interval loc= new Interval(bedrec.getContig(),bedrec.getStart(),bedrec.getEnd(),false,bedrec.getContig()+":"+bedrec.getStart()+"-"+bedrec.getEnd());
					if(dict0.getSequence(loc.getContig())==null) {
						throw new JvarkitException.ContigNotFoundInDictionary(loc.getContig(), dict0);
						}
					if(intervalTreeMap.containsOverlapping(loc)) {
						throw new IllegalArgumentException("intervals in "+this.bedPath+" overlap and it's forbiden. Intervals should have been merged. see"+loc);
						}
					intervalTreeMap.put(loc, loc);
					}
				}

			final SAMSequenceDictionary destDict = new SAMSequenceDictionary(
				intervalTreeMap.values().
					stream().
					sorted(new ContigDictComparator(dict0).createLocatableComparator()).
					map(LOC->new SAMSequenceRecord(LOC.getName(), LOC.getLengthOnReference())).
					collect(Collectors.toList())
				);
			
			
			final VCFHeader headerOut = new VCFHeader(headerIn);
			headerOut.setSequenceDictionary(destDict);
			JVarkitVersion.getInstance().addMetaData(this, headerOut);
			out.writeHeader(headerOut);
			while(iterin.hasNext()) {
				final VariantContext ctx0 = iterin.next();
				
				final Interval loc0 = intervalTreeMap.getOverlapping(ctx0).
						stream().
						filter(RGN->RGN.contains(ctx0)).
						findFirst().
						orElse(null);
				if(loc0==null) {
				    continue;
					}
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx0);
				vcb.chr(loc0.getName());
				vcb.start(1 + ctx0.getStart() - loc0.getStart());
				final int chromEnd= 1 + ctx0.getEnd() - loc0.getStart(); 
				vcb.stop(chromEnd);
				if(ctx0.hasAttribute(VCFConstants.END_KEY)) {
					vcb.attribute(VCFConstants.END_KEY,chromEnd);
					}
				out.add(vcb.make());
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}

	
	
	public static void main(final String[] args) {
		new ShiftVcf().instanceMainWithExit(args);
	}

}
