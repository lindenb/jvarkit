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
* 2015 moving to knime

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;

/**
 * 
 * 
 * 
 * @author lindenb
 *
 */

public class VcfBurdenSplitter2 
	{
	private static final Logger LOG = Logger.build().
			prefix("splitter2").
			make();
	@Parameter(names = "-ignore_filtered", description = "Ignore FILTERED variants")
	private boolean ignore_filtered=false;

	

	@Parameter(names={"-sp","--splitterName"},description="Splitter Name")
	private String splitterName = "vepso";

	@Parameter(names={"-if","--ignorefilter"},description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
	private boolean acceptFiltered = false;


	@Parameter(names={"-all_filtered","--all_filtered"},description="If defined, the group where ALL the variants are FILTERED will be saved here.")
	private File allFilteredFileOut = null;

	@Parameter(names="--maxRecordsInRam",description="Max records in RAM")
	private int maxRecordsInRam=50000;
			
	@Parameter(names="--tmpDir",description="Temporary directory")
	private File tmpDir = new File(System.getProperty("java.io.tmpdir","."));
	
	@Parameter(names="--manifestFile",description="Manifest File")
	private File manifestFile = null;

	
	/** list of available splitters */
	private final Splitter splitters[]= new Splitter[]{
			new PredictionsSplitter(),
			new SlidingWindowSplitter(1000, 500),
			new SlidingWindowSplitter(1000, 300),
			new SlidingWindowSplitter(2000, 1000),
			new SlidingWindowSplitter(2000, 500),
			new SlidingWindowSplitter(4000, 2000),
			new SlidingWindowSplitter(4000, 1000),
			new SlidingWindowSplitter(10000, 3000)/* matilde 17 Fev 2017 */
		};
	

	
	
	private MyWriter oneAndOnlyWriter=null;
	public VariantContextWriter open(final VariantContextWriter delegate) {
		if(this.oneAndOnlyWriter!=null) {
			throw new JvarkitException.ProgrammingError("This method shouldn't be invoked twice");
			}
		final Optional<Splitter> splitter=
				Arrays.asList(this.splitters).
				stream().filter(S-> splitterName.equals(S.getName())).
				findFirst();
		
		if(!splitter.isPresent()) {
			throw new JvarkitException.UserError("Cannot find a splitter named "+this.splitterName);
			}
		this.oneAndOnlyWriter = new MyWriter(delegate,splitter.get());

		return this.oneAndOnlyWriter;
		}
	
	
	private  class MyWriter extends DelegateVariantContextWriter
		{
		private final Splitter splitter;
		private String prev_contig=null;
		private SortingCollection<Interval> sortingcollection=null;
		private PrintWriter manifestWriter=null;
		private VCFInfoHeaderLine infoHeaderLine;

		MyWriter(final VariantContextWriter w,final Splitter splitter)
			{
			super(w);
			this.splitter=splitter;
			}
		
		@Override
		public void writeHeader(final VCFHeader header) {
			final VCFHeader header2= new VCFHeader(header);
			this.prev_contig=null;
			this.splitter.initialize(header);

			try {
				this.manifestWriter=new PrintWriter(VcfBurdenSplitter2.this.manifestFile);
			} catch (FileNotFoundException e) {
				throw new RuntimeIOException(e);
			}
			

			
			this.infoHeaderLine = new VCFInfoHeaderLine(VCF_HEADER_SPLITKEY, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String,"Split Names");
			header2.addMetaDataLine(this.infoHeaderLine);
			super.writeHeader(header2);
			}
		
		@Override
		public void add(final VariantContext ctx) {
			if(prev_contig==null || !ctx.getContig().equals(prev_contig)) {
				dump();
				prev_contig=ctx.getContig();
				}
			if(ctx.isFiltered() && VcfBurdenSplitter2.this.ignore_filtered)
				{
				super.add(ctx);
				return;
				}
			
			final Set<String> splitNames = getSplitNamesFor(ctx);
			if(splitNames.isEmpty())
				{
				super.add(ctx);
				return;
				}
			if(this.sortingcollection==null) {
				/* create sorting collection for new contig */
				this.sortingcollection = SortingCollection.newInstance(
						Interval.class,
						new IntervalCodec(),
						new IntervalComparator(),
						VcfBurdenSplitter2.this.maxRecordsInRam,
						VcfBurdenSplitter2.this.tmpDir
						);
				this.sortingcollection.setDestructiveIteration(true);
				}
			for(final String spltiName:splitNames)
				{
				this.sortingcollection.add(
						new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd(),false,spltiName));
				}
			
			
			super.add(new VariantContextBuilder(ctx).attribute(
					this.infoHeaderLine.getID(),
					new ArrayList<>(splitNames)).
					make());
			}
		@Override
		public void close() {
			dump();
			if(manifestWriter!=null) {
				this.manifestWriter.flush();
				if(this.manifestWriter.checkError()) {
					throw new RuntimeIOException("There was a I/O error when writing the manifest files");
				}
				this.manifestWriter.close();
			}
			super.close();
			}
		
		private Set<String> getSplitNamesFor(final VariantContext ctx){
			return this.splitter.keys(ctx);
			}

			private void dump() {
				if(this.sortingcollection==null || this.manifestWriter==null) return;
				CloseableIterator<Interval> iter;
				this.sortingcollection.doneAdding();
				iter = this.sortingcollection.iterator();
				LOG.info("dumping data for CONTIG: \""+prev_contig+"\"");
				
				final EqualRangeIterator<Interval> eqiter = new EqualRangeIterator<>(iter, new Comparator<Interval>() {
					@Override
					public int compare(final Interval o1, final Interval o2) {
						return o1.getName().compareTo(o2.getName());
						}
					});
				while(eqiter.hasNext())
					{
					final List<Interval> buffer = eqiter.next();
					
					final Interval first =  buffer.get(0);
					this.manifestWriter.print(first.getContig());
					this.manifestWriter.print('\t');
					this.manifestWriter.print(buffer.stream().map(I->I.getStart()).min((A,B)->A.compareTo(B)).get());
					this.manifestWriter.print('\t');
					this.manifestWriter.print(buffer.stream().map(I->I.getEnd()).max((A,B)->A.compareTo(B)).get());
					this.manifestWriter.print('\t');
					this.manifestWriter.print(first.getName());
					this.manifestWriter.print('\t');
					this.manifestWriter.print(buffer.size());
					this.manifestWriter.println();
					this.manifestWriter.flush();
					}
				eqiter.close();
				iter.close();iter=null;
				//dispose sorting collection
				sortingcollection.cleanup();
				sortingcollection=null;
				}

			
			}
	
	

	
	//public for knime
	public static final String VCF_HEADER_SPLITKEY="VCFBurdenSplitName";
	
	
	
	
	private static class IntervalComparator
		implements Comparator<Interval>
		{
		@Override
		public int compare(final Interval o1, final Interval o2) {
			int i = o1.getName().compareTo(o2.getName());
			if(i!=0) return i;
			
			if(!o1.getContig().equals(o2.getContig())) {
				throw new IllegalStateException("not same contig???");
			}
			i =o1.getStart() - o2.getStart();
			if(i!=0) return i;
			return o1.getEnd() - o2.getEnd();
			}
		}
	
	private static class IntervalCodec extends AbstractDataCodec<Interval>
		{
		@Override
		public Interval decode(final DataInputStream dis) throws IOException {
			String k;
			try {
				k=dis.readUTF();
			} catch(IOException err) { return null;}
			final String contig = dis.readUTF();
			final int beg= dis.readInt();
			final int end= dis.readInt();
			return new Interval(contig,beg,end,false,k);
		}
		@Override
		public void encode(final DataOutputStream dos, final Interval object) throws IOException {
			dos.writeUTF(object.getName());
			dos.writeUTF(object.getContig());
			dos.writeInt(object.getStart());
			dos.writeInt(object.getStart());
			}
		@Override
		public IntervalCodec clone() {
			return new IntervalCodec();
			}
		}

	protected boolean isDebuggingVariant(VariantContext ctx) {
		return false;
	}
	
	protected String shortName(VariantContext ctx) {
	return ctx.getContig()+":"+ctx.getStart()+":"+ctx.getAlleles();	
	}
	
	/** abstract splitter */
	private abstract class Splitter {
		public abstract Set<String> keys(final VariantContext ctx);
		public abstract String getName();
		public abstract String getDescription();
		public void initialize(final VCFHeader header) {}
		@Override
		public String toString() {
			return getName()+" "+getDescription();
			}
		}
	
	/** base vep splitter */
	private class PredictionsSplitter extends Splitter {
		private VepPredictionParser vepPredictionParser=null;
		private AnnPredictionParser annPredictionParser=null;
		PredictionsSplitter() {}
		
		
		
		@Override public String getName() { return "prediction";}
		@Override public String getDescription() { return "Variant Effect Predictions";}
		public void initialize(final VCFHeader header) {
			this.vepPredictionParser = new VepPredictionParserFactory(header).get();
			this.annPredictionParser = new AnnPredictionParserFactory(header).get();
			}
		private boolean isEmpty(final String s) {
			return s==null || s.trim().isEmpty();
		}
		@Override
		public Set<String> keys(final VariantContext ctx) {
			final Set<String> keys = new HashSet<>();
			for(final VepPrediction pred: this.vepPredictionParser.getPredictions(ctx)) {
				if(!isEmpty(pred.getFeature()))
					{
					keys.add(ctx.getContig()+"_VEP_FEATURE_"+pred.getFeature());
					}
				if(!isEmpty(pred.getGene()))
					{
					keys.add(ctx.getContig()+"_VEP_GENE_"+pred.getGene());
					}
				if(!isEmpty(pred.getSymbol()))
					{
					keys.add(ctx.getContig()+"_VEP_SYMBOL_"+pred.getSymbol());
					}
				if(!isEmpty(pred.getRefSeq()))
					{
					keys.add(ctx.getContig()+"_VEP_REFSEQ_"+pred.getRefSeq());
					}
				}
			for(final AnnPrediction pred: this.annPredictionParser.getPredictions(ctx)) {
				if(!isEmpty(pred.getFeatureType()) && !isEmpty(pred.getFeatureId())) {
					keys.add(ctx.getContig()+"_ANN_FEATURE_"+pred.getFeatureType()+"_"+pred.getFeatureId());
					}
				if(!isEmpty(pred.getGeneId()))
					{
					keys.add(ctx.getContig()+"_ANN_GENE_"+pred.getFeatureType());
					}
				}
			return keys;
			}
		}
	
	

		
	
	private class SlidingWindowSplitter extends Splitter {
		final int winsize;
		final int winshift;
		SlidingWindowSplitter(final int winsize,final int winshift)
			{
			this.winsize = winsize;
			this.winshift =winshift;
			}
		
		private int leftMostWindowStart(final VariantContext ctx) {
	    	int varstart = ctx.getStart();
	    	varstart = Math.max(0, varstart - varstart%this.winsize);
	    	while(varstart>0)
	    	{
	    		int left = varstart - this.winshift;
	    		if( left <0) left=0;
	    		if( left + this.winsize < ctx.getStart()) break;
	    		varstart = left;
	    	}
	    	return varstart;
	    }
		
		@Override
		public Set<String> keys(final VariantContext ctx)
			{
			final Set<String> set= new HashSet<>();
			int x= this.leftMostWindowStart(ctx);
			
			while(x<=ctx.getStart() && ctx.getStart()<=x+this.winsize) {
				set.add(String.format("%s_%09d_%09d", ctx.getContig(),x,x+winsize));
				x += this.winshift;
			}
	 		return set;
			}
		
		@Override
		public String getName() {
			return "split"+this.winsize+"_"+this.winshift;
			}
		@Override
		public String getDescription() {
			return "split using a sliding window of "+this.winsize+" each "+this.winshift;
			}
	}
	
	public VcfBurdenSplitter2()
		{
		}
	
	
	
	
	

	
	 
	private static class Launcher extends com.github.lindenb.jvarkit.util.jcommander.Launcher
		{
		
		@Parameter(names={"-ls","--listsplitters"},description="List available splitters and exit")
		private boolean listSplitters = false;
		@ParametersDelegate
		private VcfBurdenSplitter2 instance=new VcfBurdenSplitter2();

		
		Launcher() {
			
		
			}
		@Override
		public int doWork(final List<String> args) {
			if(this.listSplitters) {
				for(final Splitter splitter:instance.splitters) {
					stdout().println(splitter.getName()+"\t"+splitter.getDescription());
					}
				return 0;
				}
			
			if(instance.splitterName==null || instance.splitterName.isEmpty()) {
				throw new JvarkitException.CommandLineError("splitter name undefined ");
			}

			
			return doVcfToVcf(oneFileOrNull(args));
			}
		}
	
	
	public static void main(String[] args)
		{
		args=new String[]{"-h"};
		new Launcher().instanceMainWithExit(args);
		}
	}
