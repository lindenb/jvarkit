/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.vcfbed;

import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.github.lindenb.jvarkit.util.bio.bed.IndexedBedReader;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.bed.BedLine;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


/**
 * VCFBed
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
java -jar dist/vcfbed.jar -m your.bed -f '${1}|${2}|${3}|${4}&${5}'

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
	biostars=247224
	)
public class VCFBed extends Launcher
	{

	private static final Logger LOG = Logger.build(VCFBed.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-f","--format"},description="format pattern ${xx} will be replaced by column xx in the bed line. Empty lines will be ignored (no tag) but the FILTERs will be set.")
	private String formatPattern = "${1}:${2}-${3}";

	@Parameter(names={"-T","--tag"},description="use the following INFO tag name")
	private String infoName = "VCFBED";

	@Parameter(names={"-B","--bed"},description="Tribble or Tabix bed file ")
	private File tabixFile = null;

	@Parameter(names={"-m","--map"},description="unindexed bed file, will be loaded in memory (faster than tribble/tabix but memory consumming)")
	private File treeMapFile = null;

	@Parameter(names={"-fo","--filteroverlap"},description="if defined, set this as a FILTER column if one or more BED line overlap a variant")
	private String filterOverlapStr = null;

	@Parameter(names={"-fn","--filternooverlap"},description="if defined, set this as a FILTER column if not any BED line overlap a variant")
	private String filterNoOverlapStr = null;

	private IntervalTreeMap<Set<BedLine>> intervalTreeMap=null;

	
	private IndexedBedReader bedReader =null;
	private Chunk parsedFormat=null;
	
	private static abstract class Chunk
		{
		public abstract String toString(BedLine tokens);
		public Chunk next=null;
		}
	
	private static class PlainChunk extends Chunk
		{
		final String s;
		PlainChunk(final String s){this.s=s;}
		public String toString(final BedLine tokens)
			{
			return s+(next==null?"":next.toString(tokens));
			}
		@Override
		public String toString() {
			return "plain-text:\""+this.s+"\""+(this.next==null?"":";"+this.next.toString());
			}
		}
	private static class ColChunk extends Chunk
		{
		final int index;
		ColChunk(final int index){ this.index=index;}
		public String toString(final BedLine tokens)
			{
			String s= tokens.get(index);
			if(s==null) s="";
			return s+(next==null?"":next.toString(tokens));
			}
		@Override
		public String toString() {
			return "column:${"+(index+1)+"}"+(this.next==null?"":";"+this.next.toString());
			}
		}

	
	private Chunk parseFormat(final String s)
		{
		if(s==null || s.isEmpty()) return null;
		if(s.startsWith("${"))
			{
			final int j=s.indexOf('}',2);
			if(j==-1) throw new IllegalArgumentException("bad format in \""+s+"\".");
			try
				{
				final int col=Integer.parseInt(s.substring(2, j).trim());
				if(col<1) throw new IllegalArgumentException();
				final ColChunk c=new ColChunk(col-1);
				c.next=parseFormat(s.substring(j+1));
				return c;
				}
			catch(final Exception err)
				{
				 throw new IllegalArgumentException("bad format in \""+s+"\".",err);
				}
			}
		else if(s.startsWith("$"))
			{
			int j=1;
			while(j<s.length() && Character.isDigit(s.charAt(j)))
				{
				++j;
				}
			int col=Integer.parseInt(s.substring(1, j).trim());
			if(col<1) throw new IllegalArgumentException();
			ColChunk c=new ColChunk(col-1);
			c.next=parseFormat(s.substring(j));
			return c;
			}
		int i=0;
		final StringBuilder sb=new StringBuilder();
		while(i< s.length() && s.charAt(i)!='$')
			{
			sb.append(s.charAt(i));
			i++;
			}
		final PlainChunk c=new PlainChunk(sb.toString());
		c.next=parseFormat(s.substring(i));
		return c;
		}
	
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator r, VariantContextWriter w)
			 {
			try
				{
				final File srcbedfile = this.tabixFile==null?this.treeMapFile:this.tabixFile;
				final VCFHeader h2=new VCFHeader(r.getHeader());
				final VCFInfoHeaderLine infoHeader= 
						new VCFInfoHeaderLine(
								this.infoName,
								VCFHeaderLineCount.UNBOUNDED,
								VCFHeaderLineType.String,
								"metadata added from "+ srcbedfile+
								" . Format was "+this.formatPattern
								);
				
				final VCFFilterHeaderLine filterOverlap = 
						(this.filterOverlapStr==null || this.filterNoOverlapStr.trim().isEmpty()?null:
						new VCFFilterHeaderLine(this.filterOverlapStr, "Variant overlap with "+srcbedfile)	
						);
				
				final VCFFilterHeaderLine filterNoOverlap = 
						(this.filterNoOverlapStr==null || this.filterNoOverlapStr.trim().isEmpty()?null:
						new VCFFilterHeaderLine(this.filterNoOverlapStr, "Variant having NO overlap with "+srcbedfile)	
						);
				
				if(filterOverlap!=null) h2.addMetaDataLine(filterOverlap);
				if(filterNoOverlap!=null) h2.addMetaDataLine(filterNoOverlap);
				
				h2.addMetaDataLine(infoHeader);
				addMetaData(h2);
				
				final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(h2);
				w.writeHeader(h2);
				while(r.hasNext())
					{
					boolean found_overlap=false;
					final VariantContext ctx= progress.watch(r.next());
					final Set<String> annotations=new HashSet<String>();
					
					if(this.intervalTreeMap!=null) {
						for(final Set<BedLine> bedLines :this.intervalTreeMap.getOverlapping(new Interval(ctx.getContig(),ctx.getStart(),ctx.getEnd()))) {
							for(final BedLine bedLine:bedLines) {
							final String newannot=this.parsedFormat.toString(bedLine);
							found_overlap=true;
							if(!newannot.isEmpty())
								{
								annotations.add(VCFUtils.escapeInfoField(newannot));
								}
							}
						}
						
					}
					else
						{
						CloseableIterator<BedLine> iter = this.bedReader.iterator(
								ctx.getContig(),
								ctx.getStart()-1,
								ctx.getEnd()+1
								);
						while(iter.hasNext())
							{
							final BedLine bedLine = iter.next();
							
							if(!ctx.getContig().equals(bedLine.getContig())) continue;
							if(ctx.getStart() > bedLine.getEnd() ) continue;
							if(ctx.getEnd() < bedLine.getStart() ) continue;
			
							found_overlap=true;
		
							final String newannot=this.parsedFormat.toString(bedLine);
							if(!newannot.isEmpty())
								annotations.add(VCFUtils.escapeInfoField(newannot));
							}
						CloserUtil.close(iter);
						}
					
					final String filterToSet;
					if(found_overlap && filterOverlap!=null) {
						filterToSet = filterOverlap.getID();
						}
					else if(!found_overlap &&  filterNoOverlap!=null) {
						filterToSet = filterNoOverlap.getID();
						}
					else
						{
						filterToSet=null;
						}
					
					if(filterToSet==null && annotations.isEmpty())
						{
						w.add(ctx);
						continue;
						}
					final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
					if(!annotations.isEmpty()) {
						vcb.attribute(infoHeader.getID(), annotations.toArray());
						}
					if(filterToSet!=null) {
						vcb.filter(filterToSet);
						}
					else if(ctx.isNotFiltered())
						{
						vcb.passFilters();
						}
					w.add(vcb.make());
					if(w.checkError()) break;
					}
				progress.finish();
				return RETURN_OK;
				}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		}
	
	/** reads a Bed file and convert it to a IntervalTreeMap<Bedline> */
	private htsjdk.samtools.util.IntervalTreeMap<Set<com.github.lindenb.jvarkit.util.bio.bed.BedLine>> 
		readBedFileAsIntervalTreeMap(final java.io.File file) throws java.io.IOException
		{
		java.io.BufferedReader r=null;
		try
			{
			final  htsjdk.samtools.util.IntervalTreeMap<Set<com.github.lindenb.jvarkit.util.bio.bed.BedLine>> intervals = new
					 htsjdk.samtools.util.IntervalTreeMap<>();
			r=com.github.lindenb.jvarkit.io.IOUtils.openFileForBufferedReading(file);
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
	public int doWork(final List<String> args) {
		try
			{
			if(this.tabixFile==null && this.treeMapFile==null)
				{
				LOG.error("Undefined tabix or memory file");
				return -1;
				}
			else if(this.tabixFile!=null && this.treeMapFile!=null)
				{
				LOG.error("You cannot use both options: tabix/in memory bed");
				return -1;
				}
			else if( this.tabixFile!=null) {
				LOG.info("opening Bed "+this.tabixFile);
				this.bedReader= new IndexedBedReader(this.tabixFile);
				}
			else 
				{
				try {
					this.intervalTreeMap = this.readBedFileAsIntervalTreeMap(this.treeMapFile);
					LOG.info("Number of items in "+this.treeMapFile+" "+this.intervalTreeMap.size());
				}
				catch(final Exception err) {
					LOG.error(err);
					return -1;
					}
				}
			
			if(this.infoName==null || this.infoName.trim().isEmpty())
				{
				LOG.error("Undefined INFO name.");
				return -1;
				}
			
			LOG.info("parsing "+this.formatPattern);
			this.parsedFormat=parseFormat(formatPattern);
			if(this.parsedFormat==null) this.parsedFormat=new PlainChunk("");
			LOG.info("format for "+this.formatPattern+" :"+this.parsedFormat);
			return doVcfToVcf(args, outputFile);
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.bedReader);
			this.bedReader = null;
			this.intervalTreeMap=null;
			this.parsedFormat = null;			
			}
		}

	
	
	public static void main(String[] args) throws Exception
		{
		new VCFBed().instanceMainWithExit(args);
		}
}
