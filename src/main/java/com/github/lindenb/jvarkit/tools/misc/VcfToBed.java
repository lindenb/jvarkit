/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.bio.fasta.ReferenceFileSupplier;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;

/*
 
BEGIN_DOC

## Motivation

I'm lazy about using awk or bioalcidaejdk for this task and I want something that uses INFO/CIPOS and INFO/CIEND for structural variants

## Input

input is one or more VCF file

one file ending with '.list' is interpreted as a list of paths (one per lines)

if there is no input, the program reads vcf from stdin

##Example


```
$ wget -q -O - "https://github.com/hall-lab/cshl_sv_2014/blob/master/supplemental/NA12878.lumpy.vcf?raw=true" |\
	grep -A 10 '#CHROM'
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
1	869423	1	G	<DEL>	345.50	.	SVTYPE=DEL;SVLEN=-857;END=870280;STR=+-:25;IMPRECISE;CIPOS=-1,34;CIEND=0,0;EVENT=1;SUP=25;PESUP=25;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	1/1:25:25:0:0.00:71:19:52:345.50:-36,-8,-1
1	1588585	5	A	<DUP>	0.00	.	SVTYPE=DUP;SVLEN=65356;END=1653941;STR=-+:7;IMPRECISE;CIPOS=-126,1;CIEND=-2,67;EVENT=5;SUP=7;PESUP=7;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:7:7:0:74.40:139:125:13:0.00:-1,-8,-24
1	1594964	6	C	<DUP>	0.00	.	SVTYPE=DUP;SVLEN=65855;END=1660819;STR=-+:8;IMPRECISE;CIPOS=-81,2;CIEND=-1,127;EVENT=6;SUP=8;PESUP=8;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:8:8:0:77.96:153:137:15:0.00:-1,-9,-25
1	2566176	7	A	<DEL>	121.20	.	SVTYPE=DEL;SVLEN=-418;END=2566594;STR=+-:14;IMPRECISE;CIPOS=-2,68;CIEND=0,0;EVENT=7;SUP=14;PESUP=14;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:14:14:0:0.00:78:44:33:121.20:-13,-1,-12
1	2911548	8	G	<DEL>	440.34	.	SVTYPE=DEL;SVLEN=-302;END=2911850;STR=+-:20;CIPOS=0,0;CIEND=0,0;EVENT=8;SUP=20;PESUP=8;SRSUP=12;EVTYPE=PE,SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:20:8:12:0.00:185:86:99:440.34:-48,-4,-15
1	2919034	9	G	<DEL>	289.83	.	SVTYPE=DEL;SVLEN=-332;END=2919366;STR=+-:22;CIPOS=0,0;CIEND=0,0;EVENT=9;SUP=22;PESUP=10;SRSUP=12;EVTYPE=PE,SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:22:10:12:0.00:160:86:74:289.83:-31,-2,-20
1	5447229	14	G	<DUP>	380.12	.	SVTYPE=DUP;SVLEN=210;END=5447439;STR=-+:11;CIPOS=0,0;CIEND=0,0;EVENT=14;SUP=11;PESUP=1;SRSUP=10;EVTYPE=PE,SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	1/1:11:1:10:0.00:197:104:93:380.12:-39,-7,-1
1	5876603	15	G	<DEL>	0.00	.	SVTYPE=DEL;SVLEN=-928;END=5877531;STR=+-:8;CIPOS=0,0;CIEND=0,0;EVENT=15;SUP=8;PESUP=0;SRSUP=8;EVTYPE=SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:8:0:8:294.07:169:168:0:0.00:-8,-37,-117
1	5877530	16	T	<DEL>	63.31	.	SVTYPE=DEL;SVLEN=-72;END=5877602;STR=+-:13;CIPOS=0,0;CIEND=0,0;EVENT=16;SUP=13;PESUP=0;SRSUP=13;EVTYPE=SR;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:13:0:13:0.00:188:136:51:63.31:-10,-4,-54
1	6619067	19_1	T	[1:6619506[T	0.00	.	SVTYPE=BND;STR=--:7;IMPRECISE;CIPOS=-88,1;CIEND=-26,2;MATEID=19_2;EVENT=19;SUP=7;PESUP=7;SRSUP=0;EVTYPE=PE;PRIN	GT:SUP:PE:SR:GQ:DP:RO:AO:SQ:GL	0/1:7:7:0:127.76:131:117:13:0.00:-1,-14,-66


$ wget -q -O - "https://github.com/hall-lab/cshl_sv_2014/blob/master/supplemental/NA12878.lumpy.vcf?raw=true" |\
	java -jar dist/vcf2bed.jar |\
	head

1	869421	870280	1	345
1	1588458	1654008	5	0
1	1594882	1660946	6	0
1	2566173	2566594	7	121
1	2911547	2911850	8	440
1	2919033	2919366	9	289
1	5447228	5447439	14	380
1	5876602	5877531	15	0
1	5877529	5877602	16	63
1	6618978	6619069	19_1	0
```

END_DOC

*/
@Program(name="vcf2bed",
	description="vcf to bed",
	keywords={"bed","vcf"}
)
public class VcfToBed  extends Launcher {
	private enum OutputFormat {
		bed,interval
	}
	
	private static final Logger LOG = Logger.build(VcfToBed.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-R","--reference"},description="Convert the contigs of the VCF on the fly using an indexed genome. " +
			ReferenceFileSupplier.OPT_DESCRIPTION,converter=ReferenceFileSupplier.StringConverter.class)
	private ReferenceFileSupplier referenceFileSupplier=null;
	@Parameter(names={"-c","--no-ci"},description="For structural variant, ignore the extention of the boundaries using INFO/CIPOS and INFO/CIEND")
	private boolean ignoreCi = false;
	@Parameter(names={"-x","--slop"},description="Extends interval by 'x' bases on both sides. "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
	private int slopSize = 0;
	@Parameter(names={"-F","--format"},description="output format")
	private OutputFormat outputFormat = OutputFormat.bed;
	@Parameter(names={"-m","--min"},description="Optional filter: min sequence length. " + DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class)
	private String minLengthStr = null;
	@Parameter(names={"-M","--max"},description="Optional filter: max sequence length. " + DistanceParser.OPT_DESCRIPTION,splitter=NoSplitter.class)
	private String maxLengthStr = null;
	@Parameter(names={"-header","--header"},description="Print Header")
	private boolean printHeader= false;
	
	private final Set<String> contigsNotFound = new TreeSet<>();
	private SAMSequenceDictionary samSequenceDictionary = null;
	private int minLength = -1;
	private int maxLength = -1;
	
	
	private void scan(String uriOrNull,final PrintWriter pw)
		{
		VCFIterator iter = null;
		try {
			iter =  super.openVCFIterator(uriOrNull);
			final SAMSequenceDictionary dictIn = iter.getHeader().getSequenceDictionary();
			
			final ContigNameConverter ctgNameConverter ;
			final Function<String, SAMSequenceRecord> funGetSSR;
			if(this.samSequenceDictionary!=null)
				{
				ctgNameConverter=ContigNameConverter.fromOneDictionary(this.samSequenceDictionary);
				funGetSSR = C->samSequenceDictionary.getSequence(C);
				}
			else if(dictIn!=null)
				{
				ctgNameConverter=ContigNameConverter.fromOneDictionary(dictIn);
				funGetSSR = C->dictIn.getSequence(C);
				}
			else
				{
				ctgNameConverter = null;
				funGetSSR = null;
				}
			
			
			
			while(iter.hasNext())
				{
				final VariantContext ctx = iter.next();
				
				int start = ctx.getStart();
				int end = ctx.getEnd();
				
				
				
				if(!this.ignoreCi && (ctx.hasAttribute("CIPOS") || ctx.hasAttribute("CIEND"))) {
					if(ctx.hasAttribute("CIPOS")) {
						int x=0;
						try {
							x= ctx.getAttributeAsIntList("CIPOS", 0).get(0);
							}
						catch(final Throwable err)
							{
							LOG.error(err);
							x=0;
							}
						start += x;
						}
					if(ctx.hasAttribute("CIEND")) {
						int x=0;
						try {
							x= ctx.getAttributeAsIntList("CIEND", 0).get(1);
							}
						catch(final Throwable err)
							{
							LOG.error(err);
							x=0;
							}
						end += x;
						}
					//it happens for SV=BND
					if(start>end)
						{
						final int tmp = start;
						start = end;
						end = tmp;
						}
					}
				
				if(this.slopSize>0)
					{
					start -= this.slopSize;
					end += this.slopSize;
					}
				
				final String newtContig ;
				if(ctgNameConverter!=null)
					{
					newtContig = ctgNameConverter.apply(ctx.getContig());
					if(StringUtil.isBlank(newtContig)) {
						this.contigsNotFound.add(ctx.getContig());
						continue;
						}
					}
				else
					{
					newtContig = ctx.getContig();
					}
				
				if(funGetSSR!=null) {
					final SAMSequenceRecord ssr = funGetSSR.apply(newtContig);
					if(ssr!=null)
						{
						end = Math.min(ssr.getSequenceLength(),end);
						}
					}
				
				start = Math.max(1, start);
				
				final Interval interval = new Interval(newtContig,start,end);
				
				if(interval.getStart()>interval.getEnd())
					{
					LOG.info("negative interval "+interval +" for "+ctx+" in "+ uriOrNull);
					continue;
					}
				
				if(this.minLength!=-1 && interval.getLengthOnReference()<this.minLength) continue;
				if(this.maxLength!=-1 && interval.getLengthOnReference()>this.maxLength) continue;
				
				
				
				
				switch(this.outputFormat)
					{
					case bed:
						{
						pw.print(interval.getContig());
						pw.print("\t");
						pw.print(interval.getStart()-1);
						pw.print("\t");
						pw.print(interval.getEnd());
						pw.print("\t");
						pw.print(ctx.hasID()?
								ctx.getID():
								ctx.getAlleles().stream().map(A->A.getDisplayString()).collect(Collectors.joining(","))
								);
						pw.print("\t");
						pw.print(ctx.hasLog10PError()?Math.min(1000,(int)ctx.getPhredScaledQual()):".");
						/*
						pw.print("\t");
						pw.print("+");
						pw.print("\t");
						pw.print(interval.getStart()-1);
						pw.print("\t");
						pw.print(interval.getEnd());
						pw.print("\t");
						pw.print("0,0,0");
						*/
						pw.println();
						break;
						}
					case interval:
						{
						pw.print(interval.getContig());
						pw.print(":");
						pw.print(interval.getStart());
						pw.print("-");
						pw.print(interval.getEnd());
						pw.println();
						break;
						}
					default:
						{
						throw new IllegalStateException(""+this.outputFormat);
						}	
					}
			
				
				}
			}
		catch(final Exception err)
			{
			LOG.error(err);
			}
		finally
			{
			CloserUtil.close(iter);
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		if(!StringUtil.isBlank(this.maxLengthStr)) {
			this.maxLength = new DistanceParser.StringConverter().convert(this.maxLengthStr);
			}
		if(!StringUtil.isBlank(this.minLengthStr)) {
			this.minLength = new DistanceParser.StringConverter().convert(this.minLengthStr);
			}
		
		PrintWriter pw=null;
		try {
			if(this.referenceFileSupplier!=null)
				{
				final File fai = this.referenceFileSupplier.getRequired();
				this.samSequenceDictionary = SequenceDictionaryUtils.extractRequired(fai);
				}
			
			pw = super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			
			if(printHeader) {
				switch(this.outputFormat)
					{
					case bed:
						{
						pw.println("track name=vcf2bed type=bed description=\"__DESCRIPTION__\"");
						break;
						}
					case interval:
						{
						if(this.samSequenceDictionary!=null)
							{
							final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
							codec.encode(pw, new SAMFileHeader(this.samSequenceDictionary));
							}
						break;
						}
					default:
						{
						throw new IllegalStateException(""+this.outputFormat);
						}	
					}
				}
				
			
			if(args.size()==1 && args.get(0).endsWith(".list"))
				{
				final PrintWriter finalpw = pw;
				Files.lines(Paths.get(args.get(0))).filter(L->!StringUtil.isBlank(L)).forEach(
					L->{
						scan(L,finalpw);
					});
				}
			else if(args.isEmpty())
				{
				LOG.info("reading vcf from stdin");
				scan(null,pw);
				}
			else
				{
				final PrintWriter finalpw = pw;
				args.stream().forEach(L->{
					scan(L,finalpw);
					});
				}	
			if(!this.contigsNotFound.isEmpty())
				{	
				LOG.warn("The following contigs: "+String.join(",",this.contigsNotFound)+
						" were not found in the dictionaries."
						);
				}
			pw.flush();
			pw.close();
			pw=null;
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(pw);
			}
		}
	
	public static void main(final String[] args) {
		new VcfToBed().instanceMainWithExit(args);
	}

}
