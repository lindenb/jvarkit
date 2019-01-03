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


History:
* 2017 creation

*/
package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.VariantContextUtils.JexlVCMatchExp;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

@DocumentedGATKFeature(
		summary="Annotate Variants using a sliding window.",
		groupName = HelpConstants.DOCS_CAT_VARMANIP,
		extraDocs = {CommandLineGATK.class}
		)
public class WindowVariants extends AbstractVariantProcessor
	{
    @Argument(shortName="select", doc="Optional Jexl expression to use when selecting the data", required=false)
    public ArrayList<String> selectExpressions = new ArrayList<>();
    @Argument(fullName="windowShift",shortName="shift", doc="Window shift (in bp.)", required=false)
    protected int window_shift = 50;
    @Argument(fullName="windowSize",shortName="wsize", doc="Window Size (in bp.)", required=false)
    protected int window_size = 150;
    @Argument(fullName="windowName",shortName="wname", doc="INFO Attribute name that will be added", required=false)
    protected String winName = "WINDOW";
    @Argument(fullName="noemptywin",shortName="noemptywin", doc="Don't print Windows in INFO having zero match.", required=false)
    protected boolean hideZeroHitWindow = false;
    @Argument(fullName="best",shortName="best", doc="Only print the window with the hightest number of matches", required=false)
    protected boolean bestMatch = false;

    
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
    		return this.getStart()+ WindowVariants.this.window_size;
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
    	final long id=WindowVariants.this.id_generator++;
    	final VariantContext ctx;
    	final Set<Window> windows=new TreeSet<>();
    	
    	Variant(final VariantContext ctx) {
    		this.ctx = ctx;
    		}
    	
    	VariantContext build() {
			final List<String> winStrs = new ArrayList<>(this.windows.size());
			Window best=null;
			for(final Window w:this.windows) 
				{
				if(w.count_matching==0 && hideZeroHitWindow) continue;
				if( WindowVariants.this.bestMatch ) 
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

			
			
    		if((WindowVariants.this.bestMatch && best==null) )
    			{
    			return this.ctx;
    			}
    		else if(WindowVariants.this.bestMatch) {
    			return new VariantContextBuilder(this.ctx).
    					attribute(winName, best.toString()).
    					make();
    			}
    		else if( winStrs.isEmpty()) {
    			return this.ctx;
    			}
    		else
    			{    			
    			return new VariantContextBuilder(this.ctx).
    					attribute(winName, winStrs).
    					make();
    			}
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
    private List<JexlVCMatchExp> jexls = new ArrayList<>();
    private final List<Variant> variantBuffer = new ArrayList<>();
    private final Map<Window,Window> windowMap = new HashMap<>();

    
    public WindowVariants(){
    }
    
    @Override
    public void initialize() {
    	if(this.window_size<=0) {
    		throw new UserException("Bad window size.");
    	}
    	if(this.window_shift<=0) {
    		throw new UserException("Bad window shift.");
    	}
    	if(this.winName==null || this.winName.isEmpty()) {
    		throw new UserException("Bad INFO ID windowName");
    	}
    	final Map<String,String> exprMap=new HashMap<>();
    	for(final String expStr:this.selectExpressions) {
    		exprMap.put("expr"+(1+exprMap.size()), expStr);
    		}
    	this.jexls = VariantContextUtils.initializeMatchExps(exprMap);

    	final VCFHeader header= new VCFHeader(super.getVcfHeader());
    	if(header.getInfoHeaderLine(this.winName)!=null) {
    		throw new UserException("VCF header already contains the INFO header ID="+this.winName);
    	}
    	header.addMetaDataLine(new VCFInfoHeaderLine(
    			this.winName,
    			VCFHeaderLineCount.UNBOUNDED,
    			VCFHeaderLineType.String,
    			"Window : start|end|number-of-matching-variants|number-of-non-matching-variants"
    			));
    	super.writer.writeHeader(header);
    	super.initialize();
    	}
   
    
    private void flushVariants(final VariantContext ctx)
    	{
    	if(ctx==null)
    		{
    		for(final Variant v:this.variantBuffer) this.writer.add(v.build());
    		this.variantBuffer.clear();
    		this.windowMap.clear();
    		return;
    		}
		final int win_start = leftMostWindowStart(ctx);

		/* not same chromosome : dump all */
		if( !variantBuffer.isEmpty() &&
			!variantBuffer.get(0).ctx.getContig().equals(ctx.getContig()) ) {
    		for(final Variant v:this.variantBuffer) this.writer.add(v.build());
    		this.variantBuffer.clear();
    		this.windowMap.clear();
			}
    	int i=0;
    	while(i< this.variantBuffer.size())
    		{
    		final Variant curr= variantBuffer.get(i);
    		
    		if(curr.ctx.getEnd()< win_start)
    			{
    			this.writer.add(curr.build());
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
    protected VariantContext mapVariant(final VariantContext ctx,
    		final RefMetaDataTracker tracker,
    		final ReferenceContext ref,
    		final AlignmentContext context)
    	{
    	flushVariants(ctx);

    	final Variant variant = new Variant(ctx);
    	this.variantBuffer.add(variant);
    	boolean match=true;
    	if(!this.jexls.isEmpty()) {
    		match = false;
    		for(final JexlVCMatchExp exp:this.jexls)
	    		{
		    	if(VariantContextUtils.match(ctx, exp)) {
		    		match=true;
		    		break;
		    		}
	    		}
	    	}
    	
    	
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
    		if(match)
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
    public void onTraversalDone(final Long result) {
    	flushVariants(null);
       super.onTraversalDone(result);
    }

}
