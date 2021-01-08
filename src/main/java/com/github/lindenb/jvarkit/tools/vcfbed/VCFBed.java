/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.IOException;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;

import org.apache.commons.jexl2.JexlContext;

import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.jcommander.converter.FractionConverter;
import com.github.lindenb.jvarkit.jexl.JexlToString;
import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.variant.variantcontext.Breakend;

import htsjdk.variant.vcf.VCFIterator;


/**
 * VCFBed
 * 
BEGIN_DOC

## Example

Map the NCBI biosystems to a BED file using the following script:     https://gist.github.com/6024788 

```
$ gunzip -c ~/ncbibiosystem.bed.gz | head
1	69091	70008	79501	106356	30	Signaling_by_GPCR
1	69091	70008	79501	106383	50	Olfactory_Signaling_Pathway
1	69091	70008	79501	119548	40	GPCR_downstream_signaling
1	69091	70008	79501	477114	30	Signal_Transduction
1	69091	70008	79501	498	40	Olfactory_transduction
1	69091	70008	79501	83087	60	Olfactory_transduction
1	367640	368634	26683	106356	30	Signaling_by_GPCR
1	367640	368634	26683	106383	50	Olfactory_Signaling_Pathway
1	367640	368634	26683	119548	40	GPCR_downstream_signaling
1	367640	368634	26683	477114	30	Signal_Transduction
```

Now, annotate a remote VCF with the data of NCBI biosystems.

```
curl "https://raw.github.com/arq5x/gemini/master/test/test1.snpeff.vcf" |\
 sed 's/^chr//' |\
 java -jar  dist/vcfbed.jar -B ~/ncbibiosystem.bed.gz -T NCBIBIOSYS  -f '($4|$5|$6|$7)' |\
 grep -E '(^#CHR|NCBI)'

##INFO=<ID=NCBIBIOSYS,Number=.,Type=String,Description="metadata added from /home/lindenb/ncbibiosystem.bed.gz . Format was ($4|$5|$6|$7)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1094PC0005	1094PC0009	1094PC0012	1094PC0013
1	69270	.	A	G	2694.18	.	AC=40;AF=1.000;AN=40;DP=83;Dels=0.00;EFF=SYNONYMOUS_CODING(LOW|SILENT|tcA/tcG|S60|305|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008);FS=0.000;HRun=0;HaplotypeScore=0.0000;InbreedingCoeff=-0.0598;MQ=31.06;MQ0=0;NCBIBIOSYS=(79501|119548|40|GPCR_downstream_signaling),(79501|106356|30|Signaling_by_GPCR),(79501|498|40|Olfactory_transduction),(79501|83087|60|Olfactory_transduction),(79501|477114|30|Signal_Transduction),(79501|106383|50|Olfactory_Signaling_Pathway);QD=32.86	GT:AD:DP:GQ:PL	./.	./.	1/1:0,3:3:9.03:106,9,0	1/1:0,6:6:18.05:203,18,0
1	69511	.	A	G	77777.27	.	AC=49;AF=0.875;AN=56;BaseQRankSum=0.150;DP=2816;DS;Dels=0.00;EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aca/Gca|T141A|305|OR4F5|protein_coding|CODING|ENST00000335137|exon_1_69091_70008);FS=21.286;HRun=0;HaplotypeScore=3.8956;InbreedingCoeff=0.0604;MQ=32.32;MQ0=0;MQRankSum=1.653;NCBIBIOSYS=(79501|119548|40|GPCR_downstream_signaling),(79501|106356|30|Signaling_by_GPCR),(79501|498|40|Olfactory_transduction),(79501|83087|60|Olfactory_transduction),(79501|477114|30|Signal_Transduction),(79501|106383|50|Olfactory_Signaling_Pathway);QD=27.68;ReadPosRankSum=2.261	GT:AD:DP:GQ:PL	./.	./.	0/1:2,4:6:15.70:16,0,40	0/1:2,2:4:21.59:22,0,40</h:pre>
```

Another example:

```
$ tabix -h dbsnp138_00-All.vcf.gz "19:58864565-58865165" | sed '/^[^#]/s/^/chr/' |\
java -jar dist/vcfbed.jar -m your.bed -e 'bed.get(0)+"|"+bed.get(1)+"|"+bed.get(2)+"|"+bed.get(3)+"&"+bed.get(4)'

##INFO=<ID=VCFBED,Number=.,Type=String,Description="metadata added from your.bed . Format was ${1}|${2}|${3}|${4}&${5}">
(...)
chr19   58864911    rs113760967 T   C   .   .   GNO;OTHERKG;R5;RS=113760967;RSPOS=58864911;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VP=0x050100020001000102000100;WGT=1;dbSNPBuildID=132
chr19   58865054    rs893183    T   C   .   .   CAF=[0.1299,0.8701];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;PH3;R5;RS=893183;RSPOS=58865054;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051f000100;WGT=1;dbSNPBuildID=86
chr19   58865068    rs893182    T   C   .   .   CAF=[0.1299,0.8701];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;PH3;R5;RS=893182;RSPOS=58865068;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051f000100;WGT=1;dbSNPBuildID=86
chr19   58865082    rs893181    A   T   .   .   CAF=[0.1295,0.8705];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;PH3;R5;RS=893181;RSPOS=58865082;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051f000100;WGT=1;dbSNPBuildID=86
chr19   58865091    rs893180    A   G   .   .   CAF=[0.1299,0.8701];COMMON=1;G5;GNO;HD;KGPROD;KGPhase1;KGPilot123;OTHERKG;R5;RS=893180;RSPOS=58865091;RV;SAO=0;SLO;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VLD;VP=0x05010002000115051e000100;WGT=1;dbSNPBuildID=86
chr19   58865112    rs188818621 C   T   .   .   CAF=[0.9954,0.004591];COMMON=1;KGPROD;KGPhase1;R5;RS=188818621;RSPOS=58865112;SAO=0;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VP=0x050000020001000014000100;WGT=1;dbSNPBuildID=135
chr19   58865164    rs80109863  C   T   .   .   CAF=[0.9949,0.005051];COMMON=1;GNO;KGPROD;KGPhase1;OTHERKG;R5;RS=80109863;RSPOS=58865164;SAO=0;SSR=0;VC=SNV;VCFBED=chr19|58864565|58865165|A1BG&58864865;VP=0x050000020001000116000100;WGT=1;dbSNPBuildID=132
```

END_DOC

 */
@Program(name="vcfbed",
	description="Transfer information from a BED to a VCF",
	keywords={"bed","vcf","annotation"},
	biostars=247224,
	creationDate="20180406",
	modificationDate="20200623"
	)
public class VCFBed extends OnePassVcfLauncher
	{

	private static final Logger LOG = Logger.build(VCFBed.class).make();

	@Parameter(names={"-e","--expr","--jexl","--format"},description="[20180124]A JEXL Expression returning a string " + JexlToString.OPT_WHAT_IS_JEXL +". "
			+ "The variable 'bed' is the current observed BedLine (see  https://github.com/lindenb/jvarkit/blob/7bddffca3899196e568fb5e1a479300c0038f74f/src/main/java/com/github/lindenb/jvarkit/util/bio/bed/BedLine.java ). "
			+ " The variable 'ctx' or 'variant' is the current observed variant."
			+ " The variable 'line' is the original bed line"
			)
	private String formatPattern = "bed.get(0)+\":\"+bed.get(1)+\"-\"+bed.get(2)";

	@Parameter(names={"-T","--tag"},description="Name for the INFO tag name")
	private String infoName = "VCFBED";

	@Parameter(names={"-B","--bed","-m","--map"},description="Tribble or Tabix bed file. Files must be indexed unless option --fast is selected.", required=true)
	private Path inputBedFile = null;

	@Parameter(names={"--fast","--memory"},description="Load files in memory (faster than tribble/tabix but memory consumming)")
	private boolean in_memory=false;

	@Parameter(names={"-ignoreFiltered","--ignoreFiltered"},description="[20171031] Ignore FILTERed Variants (should be faster)")
	private boolean ignoreFILTERed=false;

	@Parameter(names={"-extend","--extend","--within-distance"},description="Variant and BED must be within 'x' bp." + DistanceParser.OPT_DESCRIPTION, converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int within_distance=0;	
	@Parameter(names={"-mofv","--min-overlap-vcf-fraction"},description="Minimum overlap required as a fraction of VCF record. "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_overlap_vcf_fraction=  null;	
	@Parameter(names={"-mofb","--min-overlap-bed-fraction"},description="Minimum overlap required as a fraction of BED record. "+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_overlap_bed_fraction=  null;	
	@Parameter(names={"-mofr","--min-overlap-fraction"},description="Require that the minimum fraction be satisfied for VCF OR BED."+FractionConverter.OPT_DESC,converter=FractionConverter.class,splitter=NoSplitter.class)
	private Double min_overlap_both_fraction=  null;
	@Parameter(names={"--bnd"},description="for SVTYPE=BND do not try to extend the interval to the variant and the mate breakend.")
	private boolean disable_bnd_eval  = false;

	
	
	private IntervalTreeMap<Set<BedLine>> intervalTreeMap=null;
	private IndexedBedReader bedReader =null;
	private ContigNameConverter contigNameConverter = null;
	private Function<JexlContext, String> bedJexlToString;

	
	
	
	private static class BedJEXLContext
		implements JexlContext
		{
		final BedLine bedLine;
		final VariantContext ctx;
		BedJEXLContext(final BedLine bedLine,final VariantContext ctx) {
			this.bedLine = bedLine;
			this.ctx = ctx;
			}
		@Override
		public Object get(final String name) {
			if(name.equals("ctx")) return this.ctx;
			if(name.equals("variant")) return this.ctx;
			if(name.equals("bed")) return this.bedLine;
			if(name.equals("line")) return this.bedLine.join("\t");
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
	
	private boolean testFinerIntersection(
			final Locatable variant,
			final Locatable bed
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
	
	
	/** reads a Bed file and convert it to a IntervalTreeMap<Bedline> */
	private IntervalTreeMap<Set<BedLine>>  readBedFileAsIntervalTreeMap(final Path file) throws java.io.IOException
		{
		java.io.BufferedReader r=null;
		try
			{
			final  htsjdk.samtools.util.IntervalTreeMap<Set<com.github.lindenb.jvarkit.util.bio.bed.BedLine>> intervals = new
					 htsjdk.samtools.util.IntervalTreeMap<>();
			r=com.github.lindenb.jvarkit.io.IOUtils.openPathForBufferedReading(file);
			String line;
			final com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec codec = new com.github.lindenb.jvarkit.util.bio.bed.BedLineCodec();
			while((line=r.readLine())!=null) 
				{
				if(line.startsWith("#") ||  com.github.lindenb.jvarkit.util.bio.bed.BedLine.isBedHeader(line) ||  line.isEmpty()) continue; 
				final com.github.lindenb.jvarkit.util.bio.bed.BedLine bl = codec.decode(line);
				if(bl==null || bl.getStart()>bl.getEnd()) continue;
				final htsjdk.samtools.util.Interval interval= bl.toInterval();
				Set<BedLine> set = intervals.get(interval);
				if(set==null )
					{
					set = new HashSet<>();
					intervals.put(interval,set); 
					}
				set.add(bl);
				}
			return intervals;
			}
		finally
			{
			htsjdk.samtools.util.CloserUtil.close(r);
			}
		}

	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final  VCFIterator vcfin,
			final  VariantContextWriter out)
		{	
		final VCFHeader header = vcfin.getHeader();
		final VCFHeader h2=new VCFHeader(header);
		final VCFInfoHeaderLine infoHeader= 
				new VCFInfoHeaderLine(
						this.infoName,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"metadata added from "+ this.inputBedFile+
						" . Format was " +
						this.formatPattern.replaceAll("[\"\'\\\\]"," ")
						);
		
		final VCFInfoHeaderLine infoCountBasicOverlap = 
				new VCFInfoHeaderLine(this.infoName+"_N",1,VCFHeaderLineType.Integer, "Number of raw overlap within distance "+this.within_distance+"bp.")	
				;
		
		final VCFInfoHeaderLine infoCountFinerOverlap = 
				new VCFInfoHeaderLine(this.infoName+"_C",1,VCFHeaderLineType.Integer, "Number of finer overlap within distance "+inputBedFile+" and overlap bed:"+this.min_overlap_bed_fraction+" and overlap vcf: "+this.min_overlap_vcf_fraction)	
				;
		
		if(infoCountBasicOverlap!=null) h2.addMetaDataLine(infoCountBasicOverlap);
		if(infoCountFinerOverlap!=null) h2.addMetaDataLine(infoCountFinerOverlap);
		
		h2.addMetaDataLine(infoHeader);
		
		JVarkitVersion.getInstance().addMetaData(this, h2);
		out.writeHeader(h2);
		
		while(vcfin.hasNext())
			{
			final VariantContext ctx =  vcfin.next();
			final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
			
			vcb.rmAttribute(infoCountBasicOverlap.getID());
			vcb.rmAttribute(infoCountFinerOverlap.getID());
			vcb.rmAttribute(infoHeader.getID());
			
			
			if(this.ignoreFILTERed && ctx.isFiltered())
				{
				out.add(vcb.make());
				continue;
				}
			
			final String normalizedContig =  this.contigNameConverter.apply(ctx.getContig());
			if(StringUtil.isBlank(normalizedContig)) {
				out.add(ctx);
				continue;
				}
			
			
			final Optional<? extends Locatable> bndOpt = this.disable_bnd_eval?Optional.empty():Breakend.parseInterval(ctx);
			
			final SimpleInterval theInterval;
			
			if(bndOpt.isPresent()) {
				theInterval = new SimpleInterval(
						normalizedContig,
						bndOpt.get().getStart(),
						bndOpt.get().getEnd()
						);
				}
			else
				{
				theInterval = new SimpleInterval(
						normalizedContig,
						ctx.getStart(),
						ctx.getEnd()
						);
				}
			final SimpleInterval extendedInterval;
			if(this.within_distance<=0) {
				extendedInterval = theInterval;
				}
			else
				{
				extendedInterval= new SimpleInterval(
					theInterval.getContig(),
					Math.max(1,theInterval.getStart()+this.within_distance),
					theInterval.getEnd()+this.within_distance
					);
				}
			int count_basic_overlap = 0;
			int count_finer_overlap = 0;
			final Set<String> annotations = new LinkedHashSet<>();
			
			
			if(this.intervalTreeMap!=null) {
				for(final Set<BedLine> bedLines :this.intervalTreeMap.getOverlapping(extendedInterval)) {
					for(final BedLine bedLine:bedLines) {
						count_basic_overlap++;
						if(!testFinerIntersection(theInterval,bedLine)) continue;
						count_finer_overlap++;
						final String newannot= this.bedJexlToString.apply(new BedJEXLContext(bedLine,ctx));
						if(!StringUtil.isBlank(newannot))
							{
							annotations.add(VCFUtils.escapeInfoField(newannot));
							}
						}
				  }
				}
			else
				{
				try(CloseableIterator<BedLine> iter = this.bedReader.iterator(
							theInterval.getContig(),
							Math.max(0, theInterval.getStart()-1),
							theInterval.getEnd()+1
							)) {
					while(iter.hasNext())
						{
						final BedLine bedLine = iter.next();
						
						if(!theInterval.contigsMatch(bedLine)) continue;
						if(!extendedInterval.withinDistanceOf(bedLine, this.within_distance)) continue;
						count_basic_overlap++;
						if(!testFinerIntersection(theInterval,bedLine)) continue;
						count_finer_overlap++;
	
						final String newannot= this.bedJexlToString.apply(new BedJEXLContext(bedLine, ctx));
						if(!StringUtil.isBlank(newannot))
							annotations.add(VCFUtils.escapeInfoField(newannot));
						}
					}
				catch(final IOException ioe)
					{
					LOG.error(ioe);
					throw new RuntimeIOException(ioe);
					}
					
				}
				
			
			
			if(!annotations.isEmpty()) {
				vcb.attribute(infoHeader.getID(), annotations.toArray());
				}
			vcb.attribute(infoCountBasicOverlap.getID(),count_basic_overlap);
			vcb.attribute(infoCountFinerOverlap.getID(),count_finer_overlap);
				
			out.add(vcb.make());
				
			}
		out.close();
		return 0;
		}
	

	@Override
	protected int beforeVcf() {
		if(this.inputBedFile==null)
			{
			LOG.error("Undefined bed file file");
			return -1;
			}
		if(StringUtil.isBlank(this.infoName)) {
			LOG.error("Undefined INFO name.");
			return -1;
			}
		
		if(min_overlap_both_fraction!=null) {
			if(this.min_overlap_both_fraction <=0.0)
				{
				LOG.error("bad value for min_overlap_both_fraction");
				return -1;
				}
			this.min_overlap_bed_fraction = min_overlap_both_fraction;
			this.min_overlap_vcf_fraction = min_overlap_both_fraction;
			}
		if(this.min_overlap_bed_fraction!=null && this.min_overlap_bed_fraction <=0.0)
			{
			LOG.error("bad value for min_overlap_bed_fraction");
			return -1;
			}
		if(this.min_overlap_vcf_fraction!=null && this.min_overlap_vcf_fraction <=0.0)
			{
			LOG.error("bad value for min_overlap_vcf_fraction");
			return -1;
			}
		
		this.bedJexlToString = new JexlToString(this.formatPattern);

			
		if(!this.in_memory) {
			try 
				{
				this.bedReader = new IndexedBedReader(this.inputBedFile.toString());
				this.contigNameConverter = ContigNameConverter.fromContigSet(this.bedReader.getContigs());
				this.intervalTreeMap = null;
				}
			catch(final IOException err)
				{
				LOG.error(err);
				return -1;
				}
			}
		else 
			{
			try {
				this.bedReader = null;
				this.intervalTreeMap = this.readBedFileAsIntervalTreeMap(this.inputBedFile);
				this.contigNameConverter = ContigNameConverter.fromIntervalTreeMap(this.intervalTreeMap);
				}
			catch(final Exception err) {
				LOG.error(err);
				return -1;
				}
			}
		return 0;
		}
	
	@Override
	protected void afterVcf() {
		CloserUtil.close(this.bedReader);
		this.bedReader = null;
		this.intervalTreeMap=null;
		}

	
	
	public static void main(final String[] args)
		{
		new VCFBed().instanceMainWithExit(args);
		}
}
