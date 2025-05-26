/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
package com.github.lindenb.jvarkit.tools.vcfbigwig;

import java.io.Closeable;
import java.io.IOException;
import java.io.Reader;
import java.net.URL;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.jexl2.JexlContext;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bbfile.BigBedFeatureAsList;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.jexl.JexlToString;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/*
BEGIN_DOC

## Example

```bash
$ java -jar dist/vcfbigbed.jar \
	-B "http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/JASPAR2022_hg19.bb" \
	--format  'bed.get(1)+"-"+bed.get(2)+":"+bed.get(6)'  \
	input.vcf

(...)
##INFO=<ID=JASPAR2022_hg19,Number=.,Type=String,Description="Values from bigbed file: http://expdata.cmmt.ubc.ca/JASPAR/downloads/UCSC_tracks/2022/JASPAR2022_hg19.bb format bed.get(1)+\"-\"+bed.get(2)+\":\"+bed.get(6)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr22	41697508	.	A	C	48.67	.	AC=2;AN=10;BQB=0.572843;DP=36;DP4=19,7,3,5;HOB=0.32;ICB=0.425;JASPAR2022_hg19=41697498-41697509:ZBTB12,41697506-41697515:NKX2-8,41697507-41697524:Spi1,41697498-41697509:Stat5a::Stat5b,41697503-41697509:Foxn1,41697507-41697519:ZNF263,41697493-41697509:BCL6,41697503-41697517:TFAP2C,41697502-41697510:NR2C2;MQ=60;MQ0F=0;MQB=1;MQSB=1;RPB=0.658863;SGB=10.3229;VDB=0.693968	GT:PL	0/0:0,9,47	0/0:0,18,73
```


END_DOC
*/
@Program(name="vcfbigbed",
	description="Annotate a VCF with values from a bigbed file",
	keywords={"vcf","wig","wiggle","bigbed","bed"},
	creationDate="20220107",
	modificationDate="20220107",
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfBigBed extends OnePassVcfLauncher {
	private static final Logger LOG = Logger.of(VcfBigBed.class);

	
	/** describe a BigWig Resource */
	private  class BigBedResource
		implements Closeable
		{
		private final String tag;
		private final String biwWigFile;
		private final BBFileReader bbFileReader;
		private final ContigNameConverter contigNameConverter;
		private Locatable lastInterval = null;
		private final List<BigBedFeatureAsList> buffer = new ArrayList<>();

		BigBedResource(final String biwWigFile) throws IOException {
			this.biwWigFile = biwWigFile;
			
			String tag;
			LOG.info("opening "+this.biwWigFile+". Please wait...");
			if(IOUtil.isUrl(this.biwWigFile)) {
				final URL uri = new URL(this.biwWigFile);
				final String path = uri.getPath();
				tag = path.substring(path.lastIndexOf('/') + 1);
				int i= tag.indexOf('.');
				if(i>0)  tag= tag.substring(0,i);
				//final SeekableStream seek = SeekableStreamFactory.getInstance().getStreamFor(this.biwWigFile);
				//final SeekableStream buffered = SeekableStreamFactory.getInstance().getBufferedStream(seek,1_000_000);
				this.bbFileReader = new BBFileReader(this.biwWigFile);
				}
			else
				{
				final Path path = Paths.get(this.biwWigFile);
				IOUtil.assertFileIsReadable(path);
				tag = IOUtils.getFilenameWithoutCommonSuffixes(path);
				this.bbFileReader = new BBFileReader(this.biwWigFile);
				}
			this.tag = tag;
			if(StringUtil.isBlank(this.tag)) throw new JvarkitException.UserError("Bad TAG for "+this.biwWigFile);
			this.contigNameConverter = ContigNameConverter.fromContigSet(new HashSet<>(this.bbFileReader.getChromosomeNames()));
			}
		
	
		public String getTag() {
			return this.tag;
			}
		
		
		public Iterator<BigBedFeatureAsList> query(final VariantContext ctx) {
			final String variantChrom=  this.contigNameConverter.apply(ctx.getContig());
			if(StringUtils.isBlank(variantChrom)) {
				return Collections.emptyIterator();
				}
			if(this.lastInterval==null || !this.lastInterval.contains(ctx)) {
				this.buffer.clear();
				this.lastInterval = new SimpleInterval(variantChrom, ctx.getStart(), Math.max(ctx.getEnd(), ctx.getStart()+VcfBigBed.this.bigbedBufferSize));
				final Iterator<BedFeature> iter = this.bbFileReader.getBigBedIterator(
						variantChrom,
						this.lastInterval.getStart()-1,
						this.lastInterval.getContig(),
						this.lastInterval.getEnd(),
						false
						);
			
				while(iter.hasNext()) {
					this.buffer.add(new BigBedFeatureAsList(iter.next()));
					}
				}
			return this.buffer.stream().
					filter(B->B.overlaps(ctx)).
					filter(B->testFinerIntersection(ctx,B)).
					iterator();
			}

		
		
		@Override
		public void close() {
			CloserUtil.close(this.bbFileReader);
			}
		}
	
	private static class BedJEXLContext implements JexlContext {
		final BigBedFeatureAsList bedLine;
		final VariantContext ctx;
		BedJEXLContext(final BigBedFeatureAsList bedLine,final VariantContext ctx) {
			this.bedLine = bedLine;
			this.ctx = ctx;
			}
		@Override
		public Object get(final String name) {
			if(name.equals("ctx")) return this.ctx;
			if(name.equals("variant")) return this.ctx;
			if(name.equals("bed")) return this.bedLine;
			if(name.equals("line")) return String.join("\t",this.bedLine);
			return false;
			}
		@Override
		public boolean has(final String key) {
			if(key.equals("ctx")) return true;
			if(key.equals("variant")) return true;
			if(key.equals("bed")) return true;
			if(key.equals("line")) return true;
			return false;
			}
		@Override
		public void set(final String key, Object arg1) {
			throw new UnsupportedOperationException();
			}
		@Override
		public String toString() {
			return "JexlContext for BedLine "+this.bedLine;
			}
		}
	
	@Parameter(names={"-B","--bed","--bigbed"},description= "Path to the bigbed file.",required=true)
	private String userBigBedFileUri = null;
	@Parameter(names={"-T","--tag","-tag"},description="Name of the INFO tag. default: name of the bigbed")
	private String userVcfTag = null;
	@Parameter(names={"-mofv","--min-overlap-vcf-fraction"},description="Minimum overlap required as a fraction of VCF record. "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_overlap_vcf_fraction=  null;	
	@Parameter(names={"-mofb","--min-overlap-bed-fraction"},description="Minimum overlap required as a fraction of BIGBED record. "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_overlap_bed_fraction=  null;	
	@Parameter(names={"-e","--expr","--jexl","--format"},description="A JEXL Expression returning a string " + JexlToString.OPT_WHAT_IS_JEXL +". "
			+ "The variable 'bed' is the current observed Bed Line. It implements java.util.List<String> and Locatable."
			+ " The variable 'ctx' or 'variant' is the current observed variant."
			+ " The variable 'line' is the original bed line"
			)
	private String formatPattern = "bed.get(0)+\":\"+bed.get(1)+\"-\"+bed.get(2)";
	@Parameter(names={"--bufferSize"},description= "When we're looking for bed in a lare bigbed file,"
			+ " load the bed items in an interval of 'N' bases instead of doing a random access for each variant. "+
			DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.jcommander.NoSplitter.class)
	private int bigbedBufferSize= 10_000;


	private BigBedResource bigBedResource = null;
	private Function<JexlContext, String> bedJexlToString = null;
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	
	@Override
	protected int beforeVcf() {
		Reader fr =null;
			try
				{
				if(StringUtil.isBlank(this.userBigBedFileUri))
					{
					LOG.info("Undefined BgBed file ");
					return -1;
					}
				try
					{
					this.bigBedResource = new BigBedResource(this.userBigBedFileUri);
					}
				catch(final IOException err) {
					LOG.error(err);
					return -1;
					}
				
				
				
				this.bedJexlToString = new JexlToString(this.formatPattern);

				
				return 0;
				}
			catch(final Throwable err)
				{
				LOG.error(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(fr);
				}
			}
	
	private boolean testFinerIntersection(
			final VariantContext variant,
			final BigBedFeatureAsList bed
			) 
		{
		final int overlap_len =  CoordMath.getOverlap(variant.getStart(), variant.getEnd(), bed.getStart(), bed.getEnd());
		if(this.min_overlap_bed_fraction!=null)
			{
			final double bedL = bed.getLengthOnReference();
			if(bedL==0.0) return false; 
			if(overlap_len/bedL < this.min_overlap_bed_fraction) return false;
			}
		if(this.min_overlap_vcf_fraction!=null)
			{
			final double variantL = variant.getLengthOnReference();
			if(variantL==0.0) return false; 
			if(overlap_len/variantL < this.min_overlap_vcf_fraction) return false;
			}
		return true;
		}
	
	@Override
	protected int doVcfToVcf(final String inputName, final VCFIterator r, final VariantContextWriter w) {
		LOG.debug("X1");
		final VCFHeader header = r.getHeader();
		final VCFHeader h2=new VCFHeader(header);
		
		if(h2.getInfoHeaderLine(this.bigBedResource.getTag())!=null)
			{
			throw new JvarkitException.DuplicateVcfHeaderInfo(h2,this.bigBedResource.getTag());
			}
			
		h2.addMetaDataLine(new VCFInfoHeaderLine(
				this.bigBedResource.getTag(),
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String,
				"Values from bigbed file: "+this.userBigBedFileUri+ " format "+this.formatPattern
				));
			
			
		JVarkitVersion.getInstance().addMetaData(this, h2);
		w.writeHeader(h2);
		
		final Set<String> values = new HashSet<>();
		while(r.hasNext())
			{
			
			final VariantContext ctx = r.next();
			values.clear();
			
			for(final Iterator<BigBedFeatureAsList> iter=this.bigBedResource.query(ctx) ;  iter.hasNext();) {
				
				final BigBedFeatureAsList item= iter.next();
				final String s = this.bedJexlToString.apply(new BedJEXLContext(item,ctx));
				if(StringUtils.isBlank(s)) continue;
				values.add(s);
				}
			if(values.isEmpty())
				{
				w.add(ctx);
				}
			else
				{
				w.add(new VariantContextBuilder(ctx).
					attribute(this.bigBedResource.getTag(), new ArrayList<>(values)).
					make()
					);
				}
			}
		return 0;
		}
	
	@Override
	protected void afterVcf() {
		CloserUtil.close(this.bigBedResource);
		}
	
	public static void main(final String[] args) throws IOException
		{
		new VcfBigBed().instanceMainWithExit(args);
		}
}
