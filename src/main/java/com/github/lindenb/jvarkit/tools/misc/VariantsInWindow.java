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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.bio.DistanceParser;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.JexlVariantPredicate;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**

BEGIN_DOC

## Example

```
##fileformat=VCFv4.2
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=TOO_MANY_CLOSE_VARIANTS,Description="Filter defined in vcfwindowvariants">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)">
##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AF2,Number=1,Type=Float,Description="Max-likelihood estimate of the first and second group ALT allele frequency (assuming HWE)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
[lindenb@kaamelot-master01 jvarkit-git]$ java -jar dist/variantsinwindow.jar ~/src/gatk-ui/testdata/mutations.vcf --treshold 1 -shift 1 -windowSize 10
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4
rotavirus	51	.	A	G	22.55	.	AC1=2;AF1=0.25;BQB=1;DP=944;DP4=849,0,93,0;FQ=23.7972;G3=0.75,0,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.993129;SGB=-61.9012;VDB=3.53678e-05;WINDOW=41|50|1|0,42|51|1|0,43|52|1|0,44|53|1|0,45|54|1|0,46|55|1|0,47|56|1|0,48|57|1|0,49|58|1|0,50|59|1|0,51|60|1|0	GT:PL	0/0:0,255,134	0/0:0,255,127	0/0:0,255,137	1/1:70,255,0
rotavirus	91	.	A	T	5.45	.	AC1=1;AF1=0.124963;BQB=0.951201;DP=1359;DP4=1134,0,225,0;FQ=5.8713;MQ=60;MQ0F=0;MQB=1;PV4=1,4.80825e-05,1,1;RPB=0.0393173;SGB=-369.163;VDB=0.313337;WINDOW=81|90|1|0,82|91|1|0,83|92|1|0,84|93|1|0,85|94|1|0,86|95|1|0,87|96|1|0,88|97|1|0,89|98|1|0,90|99|1|0,91|100|1|0	GT:PL	0/0:0,255,133	0/1:40,0,31	0/0:0,255,134	0/0:0,255,82
rotavirus	130	.	T	C	4.12	.	AC1=1;AF1=0.124933;BQB=1;DP=1349;DP4=1139,0,204,0;FQ=4.48321;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.762964;SGB=-335.275;VDB=0.00084636;WINDOW=120|129|1|0,121|130|1|0,122|131|1|0,123|132|1|0,124|133|1|0,125|134|1|0,126|135|1|0,127|136|1|0,128|137|1|0,129|138|1|0,130|139|1|0	GT:PL	0/1:38,0,35	0/0:0,255,132	0/0:0,255,132	0/0:0,255,79
rotavirus	232	.	T	A	5.45	.	AC1=1;AF1=0.124959;BQB=1;DP=1308;DP4=1098,0,207,0;FQ=5.87117;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.453119;SGB=-340.116;VDB=0.0100544;WINDOW=222|231|1|0,223|232|1|0,224|233|1|0,225|234|1|0,226|235|1|0,227|236|1|0,228|237|1|0,229|238|1|0,230|239|1|0,231|240|1|0,232|241|1|0	GT:PL	0/1:40,0,35	0/0:0,255,135	0/0:0,255,132	0/0:0,255,81
rotavirus	267	.	C	G	4.76	.	AC1=1;AF1=0.124951;BQB=0.953186;DP=1393;DP4=1156,0,234,0;FQ=5.15288;MQ=60;MQ0F=0;MQB=1;PV4=1,5.65123e-05,1,1;RPB=0.0284076;SGB=-383.686;VDB=0.367507;WINDOW=257|266|1|0,258|267|1|0,259|268|1|0,260|269|1|0,261|270|1|0,262|271|1|0,263|272|1|0,264|273|1|0,265|274|1|0,266|275|1|0,267|276|1|0	GT:PL	0/1:39,0,32	0/0:0,255,132	0/0:0,255,136	0/0:0,255,76
rotavirus	424	.	A	G	52.99	.	AC1=1;AF1=0.125;BQB=0.956333;DP=1555;DP4=1096,206,198,55;FQ=53.5713;MQ=60;MQ0F=0;MQB=1;MQSB=1;PV4=0.0270045,4.0796e-05,1,1;RPB=0.0600948;SGB=-160.094;VDB=0.000623759;WINDOW=414|423|1|0,415|424|1|0,416|425|1|0,417|426|1|0,418|427|1|0,419|428|1|0,420|429|1|0,421|430|1|0,422|431|1|0,423|432|1|0,424|433|1|0	GT:PL	0/0:0,255,200	0/1:89,0,51	0/0:0,255,189	0/0:0,255,153
rotavirus	520	.	T	A	53.99	.	AC1=1;AF1=0.125;BQB=0.215002;DP=2372;DP4=1055,856,223,231;FQ=54.5713;MQ=60;MQ0F=0;MQB=1;MQSB=1;PV4=0.0211343,5.03431e-23,1,1;RPB=0.805995;SGB=-738.702;VDB=0.730588;WINDOW=510|519|1|0,511|520|1|0,512|521|1|0,513|522|1|0,514|523|1|0,515|524|1|0,516|525|1|0,517|526|1|0,518|527|1|0,519|528|1|0,520|529|1|0GT:PL	0/0:0,255,231	0/0:0,255,225	0/1:90,0,112	0/0:0,255,204
(...)
rotavirus	1054	.	C	G	15.65	TOO_MANY_CLOSE_VARIANTS	AC1=2;AF1=0.249999;BQB=1;DP=487;DP4=0,364,0,120;FQ=16.8692;G3=0.75,2.21169e-28,0.25;HWE=0.0339211;MQ=60;MQ0F=0;MQB=1;PV4=1,1,1,1;RPB=0.95941;SGB=42.7815;VDB=1.4013e-45;WINDOW=1044|1053|3|0,1045|1054|2|0,1046|1055|1|0,1047|1056|1|0,1048|1057|1|0,1049|1058|1|0,1050|1059|1|0,1051|1060|1|0,1052|1061|1|0,1053|1062|1|0,1054|1063|2|0	GT:PL	0/0:0,255,90	1/1:63,235,0	0/0:0,255,99	0/0:0,132,66
rotavirus	1064	.	G	A	21.56	TOO_MANY_CLOSE_VARIANTS	AC1=2;AF1=0.25;BQB=0.683886;DP=250;DP4=0,219,0,31;FQ=22.8019;G3=0.75,2.37734e-17,0.25;HWE=0.033921;MQ=60;MQ0F=0;MQB=1;PV4=1,1.22605e-06,1,1;RPB=0.935144;SGB=8.40135;VDB=2.70971e-16;WINDOW=1054|1063|2|0,1055|1064|1|0,1056|1065|1|0,1057|1066|1|0,1058|1067|1|0,1059|1068|1|0,1060|1069|1|0,1061|1070|1|0,1062|1071|1|0,1063|1072|1|0,1064|1073|1|0	GT:PL	0/0:0,244,70	0/0:0,199,65	0/0:0,217,68	1/1:69,84,0
```

END_DOC

*/
@Program(
		name="variantsinwindow",
		description="Annotate Number of Variants overlaping a sliding window.",
		keywords={"vcf","annotation"},
		biostars=291144
		)
public class VariantsInWindow extends Launcher{
	private static final Logger LOG = Logger.build(VariantsInWindow.class).make();

	 @Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	 private File outputFile = null;

	@Parameter(names={"-vf","--variant-filter"},description="Variants we want to keep. Variant FAILING that Jexl expression will be excluded from the window." +JexlVariantPredicate.PARAMETER_DESCRIPTION,converter=JexlVariantPredicate.Converter.class)
	private Predicate<VariantContext> variantFilter = JexlVariantPredicate.create("");
	@Parameter(names={"-S","--shift","--windowShift"},description="Window shift."+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
    protected int window_shift = 50;
    @Parameter(names= {"-W","--windowSize"},description="Window Size." + DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=com.github.lindenb.jvarkit.util.jcommander.NoSplitter.class)
    protected int window_size = 150;
    @Parameter(names="-windowName",description="INFO Attribute name that will be added")
    protected String winName = "WINDOW";
    @Parameter(names="-noemptywin",description="Don't print Windows in INFO having zero match.")
    protected boolean hideZeroHitWindow = false;
    @Parameter(names= {"--best","-best"}, description="Only print the window with the hightest number of matches")
    protected boolean bestMatch = false;
    @Parameter(names={"-filter","--filter"}, description="if --treshold is != -1 and the number of matches is greater than threshold, set this FILTER")
    protected String filterName = "TOO_MANY_CLOSE_VARIANTS";
    @Parameter(names={"-treshold","--treshold"}, description="Number of variants to set the FILTER")
    protected int treshold = -1;

    
    /** a sliding window */
    private class Window
    	implements Locatable,Comparable<Window>
    	{
    	final String contig;
    	final int start;
    	int count_matching=0;
    	int count_notmatching=0;
    	Window(final String contig, final int start) {
    		this.contig = contig;
    		this.start=start;
    		}
    	@Override
    	public String getContig() {
    		return contig;
    		}
    	@Override
    	public int getStart() {
    		return this.start;
    		}
    	@Override
    	public int getEnd() {
    		return this.getStart()+ VariantsInWindow.this.window_size -1;
    		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((contig == null) ? 0 : contig.hashCode());
			result = prime * result + start;
			return result;
		}
		
		@Override
		public boolean equals(final Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			final Window other = (Window) obj;
			if (start != other.start) {
				return false;
			}
			if (contig == null) {
				if (other.contig != null) {
					return false;
				}
			} else if (!contig.equals(other.contig)) {
				return false;
			}
			return true;
		}
		
		public int compareTo(final Window w) 
			{
			int i= this.contig.compareTo(w.contig);
			if(i!=0) return i;
			return this.start - w.start;
			}
		
		@Override
		public String toString() {
			return new StringBuilder().
					append(Math.max(0, getStart())).append("|").
					append(getEnd()).append("|").
					append(count_matching).append("|").
					append(count_notmatching).
					toString();
			}
    	
    	}
    
    /** wrapper arount a variant context , keep track of the associated windows */
    private class Variant
    	{
    	final long id=VariantsInWindow.this.id_generator++;
    	final VariantContext ctx;
    	final Set<Window> windows=new TreeSet<>();
    	
    	Variant(final VariantContext ctx) {
    		this.ctx = ctx;
    		}
    	/** build a new variant to output */
    	VariantContext build() {
    		boolean set_filter=false;
			final List<String> winStrs = new ArrayList<>(this.windows.size());
			Window best=null;
			for(final Window w:this.windows) 
				{
				if(w.count_matching==0 && hideZeroHitWindow) continue;
				
				
				if( VariantsInWindow.this.treshold>=0 && w.count_matching > VariantsInWindow.this.treshold ) {
    				set_filter = true;
    				}
				
				if( VariantsInWindow.this.bestMatch ) 
					{
					if(best==null || best.count_matching < w.count_matching)
						{
						best= w;
						}
					
					}
				else
					{
					winStrs.add(w.toString());
					}
				}
			final VariantContextBuilder vcb;
			
    		if((VariantsInWindow.this.bestMatch && best==null) )
    			{
    			return this.ctx;
    			}
    		else if(VariantsInWindow.this.bestMatch) {
    			vcb = new VariantContextBuilder(this.ctx).
    					attribute(winName, best.toString())
    					;
    			
    			}
    		else if( winStrs.isEmpty()) {
    			return this.ctx;
    			}
    		else
    			{    			
    			vcb= new VariantContextBuilder(this.ctx).
    					attribute(winName, winStrs);
    			}
    		if(set_filter ) 
    			{
				vcb.filter(VariantsInWindow.this.filterName);
				}
    		
    		return vcb.make();
    		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + (int) (id ^ (id >>> 32));
			return result;
		}
		
		@Override
		public boolean equals(final Object obj) {
			if (this == obj) {
				return true;
			}
			if (obj == null) {
				return false;
			}
			return (id == Variant.class.cast(obj).id);
			}
    	
    	}
    
    private long id_generator=0L;
    private final List<Variant> variantBuffer = new ArrayList<>();
    private final Map<Window,Window> windowMap = new HashMap<>();


    
    private void flushVariants(final VariantContextWriter writer,final VariantContext ctx)
    	{
    	if(ctx==null)
    		{
    		this.variantBuffer.stream().forEach(V->writer.add(V.build()));
    		this.variantBuffer.clear();
    		this.windowMap.clear();
    		return;
    		}
		final int win_start = leftMostWindowStart(ctx);

		/* not same chromosome : dump all */
		if( !variantBuffer.isEmpty() &&
			!variantBuffer.get(0).ctx.getContig().equals(ctx.getContig()) ) {
			this.variantBuffer.stream().forEach(V->writer.add(V.build()));
    		this.variantBuffer.clear();
    		this.windowMap.clear();
			}
    	int i=0;
    	while(i< this.variantBuffer.size())
    		{
    		final Variant curr= variantBuffer.get(i);
    		
    		if(curr.ctx.getEnd()< win_start)
    			{
    			writer.add(curr.build());
    			this.variantBuffer.remove(i);
    			
    			for(final Window win:curr.windows)
	    			{
	    			if(win.getEnd()< win_start)
	    				{
	    				this.windowMap.remove(win);
	    				}
	    			}
    			}
    		else
    			{
    			++i;
    			}
    		}
    	}
    
    private int leftMostWindowStart(final VariantContext ctx) {
    	int varstart = ctx.getStart();
    	varstart = varstart - varstart%this.window_size;
    	while(varstart>0)
    	{
    		int left = varstart - this.window_shift;
    		if( left + this.window_size < ctx.getStart()) break;
    		varstart = left;
    	}
    	return varstart;
    }
    
    /** this method is called my mapToMany */
    private VariantContext mapVariant(final VariantContextWriter writer,final VariantContext ctx)
    	{
    	flushVariants(writer,ctx);

    	final Variant variant = new Variant(ctx);
    	this.variantBuffer.add(variant);
    	final boolean rejectVariant = !this.variantFilter.test(ctx);
    	
    	
    	
    	int chromStart=  leftMostWindowStart(ctx);
    	while(!(chromStart+this.window_size < ctx.getStart() ||
    			chromStart > ctx.getEnd()
    			))
    		{
    		Window win = new Window(ctx.getContig(), chromStart);
    		
    		if(this.windowMap.containsKey(win)) {
    			win = this.windowMap.get(win);
    		} else
    			{
    			 this.windowMap.put(win,win);
    			}
    		if(!rejectVariant)
    			{
    			win.count_matching++;
    			}
    		else
    			{
    			win.count_notmatching++;
    			}
    		variant.windows.add(win);
    		chromStart+=this.window_shift;
    		}
	    	
    	return null;
    	}

    
    @Override
    protected int doVcfToVcf(
    		final String inputName,
    		final VCFIterator in,
    		final VariantContextWriter writer
    		) {
    	final VCFHeader header= new VCFHeader(in.getHeader());
    	if(header.getInfoHeaderLine(this.winName)!=null) {
    		LOG.error("VCF header already contains the INFO header ID="+this.winName);
    		}

    	header.addMetaDataLine(
			this.bestMatch?
    			new VCFInfoHeaderLine(
	    			this.winName,
	    			1,
	    			VCFHeaderLineType.String,
	    			"Window : start|end|pass-variants|filter-variants"
	    			):
	    		new VCFInfoHeaderLine(
	    			this.winName,
	    			VCFHeaderLineCount.UNBOUNDED,
	    			VCFHeaderLineType.String,
	    			"Window : start|end|pass-variants|filter-variants"
	    			));
    	
    	if(!StringUtil.isBlank(this.filterName) && this.treshold>0)
    		{
    		if(header.getInfoHeaderLine(this.filterName)!=null) {
        		LOG.error("VCF header already contains the FORMAT header ID="+this.filterName);
        		}
    		header.addMetaDataLine(new VCFFilterHeaderLine(
        			this.filterName,
        			"Filter defined in "+getProgramName()
        			));
    		}
    	final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().logger(LOG).dictionary(header).build();
    	writer.writeHeader(header);
    	while(in.hasNext()) {
    		final VariantContext ctx = progress.apply(in.next());
    		mapVariant(writer,ctx);
    		}
    	flushVariants(writer,null);
    	progress.close();
    	return 0;
    	}
    
    @Override
    public int doWork(final List<String> args) {
    	if(this.window_size<=0) {
    		LOG.error("Bad window size.");
    		return -1;
    	}
    	if(this.window_shift<=0) {
    		LOG.error("Bad window shift.");
    		return -1;
    	}
    	if(StringUtil.isBlank(this.winName)) {
    		LOG.error("Bad INFO ID windowName");
    		return -1;
    		}
        return doVcfToVcf(args, this.outputFile);
    	}
    	
    	
    public static void main(final String[] args) {
		new VariantsInWindow().instanceMainWithExit(args);
	}

}
