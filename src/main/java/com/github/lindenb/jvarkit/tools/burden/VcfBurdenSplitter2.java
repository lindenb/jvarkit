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
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

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

	@Parameter(names={"-vepFeature","--vepFeature"},description="enable VEP 'FEATURE' (transcript)")
	private boolean enableVepFeature = false;

	@Parameter(names={"-vepHgnc","--vepHgnc"},description="enable VEP 'HGNC'")
	private boolean enableVepHgnc = false;

	@Parameter(names={"-vepEnsg","--vepEnsg"},description="enable VEP 'ENSG'")
	private boolean enableVepEnsg = false;

	@Parameter(names={"-vepEnst","--vepEnst"},description="enable VEP 'FEATURE' starting with 'ENST'")
	private boolean enableVepEnst = false;

	@Parameter(names={"-vepEnsp","--vepEnsp"},description="enable VEP 'ENSP'")
	private boolean enableVepEnsp = false;

	@Parameter(names={"-vepSymbol","--vepSymbol"},description="enable VEP 'SYMBOL'")
	private boolean enableVepSymbol = false;

	@Parameter(names={"-vepRefSeq","--vepRefSeq"},description="enable VEP 'SYMBOL'= XM_ or NM_")
	private boolean enableVepRefSeq = false;

	@Parameter(names={"-all_nm","--all_nm"},description="enable grouping by ALL_NM : gene not empty and transcript starting with NM_")
	private boolean enableAllNM = false;

	@Parameter(names={"-all_refseq","--all_refseq"},description="enable grouping by ALL_REFSEQ: gene not empty and transcript NOT starting with ENST")
	private boolean enableAllRefSeq = false;

	@Parameter(names={"-all_enst","--all_enst"},description="enable grouping by ALL_ENST: gene starting with ENST")
	private boolean enableAllEnst = false;

	@Parameter(names={"-all_transcripts","--all_transcripts"},description="enable grouping by all transcript for a gene using gene name (e.g Nm_12345)")
	private boolean enableAllTranscript = false;

	@Parameter(names={"-all_genes","--all_genes"},description="enable grouping by all transcript for a gene using transcript name (e.g PRKCB1)")
	private boolean enableAllGenes = false;

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
			new VepSplitter(),
			new SoVepSplitter(),
			new ZeroVepFilter(),
			new HighDamageVepSplitter(),
			new SlidingWindowSplitter(1000, 500),
			new SlidingWindowSplitter(1000, 300),
			new SlidingWindowSplitter(2000, 1000),
			new SlidingWindowSplitter(2000, 500),
			new SlidingWindowSplitter(4000, 2000),
			new SlidingWindowSplitter(4000, 1000),
			new SlidingWindowSplitter(10000, 3000)/* matilde 17 Fev 2017 */
		};
	

	
	
	private MyWriter oneAndOnlyWriter=null;
	public VariantContextWriter open(VariantContextWriter delegate) {
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
	private class VepSplitter extends Splitter {
		private VepPredictionParser vepPredictionParser=null;
		VepSplitter() {}
		
		
		
		@Override public String getName() { return "vep";}
		@Override public String getDescription() { return "Ensembl Variant Effect Prediction";}
		public void initialize(final VCFHeader header) {
			this.vepPredictionParser = new VepPredictionParserFactory(header).get();
			}
		public boolean accept(final VepPrediction pred,final VariantContext origin) {
			return true;
		}
		private boolean isEmpty(final String s) {
			return s==null || s.trim().isEmpty();
		}
		@Override
		public Set<String> keys(final VariantContext ctx) {
			final Set<String> keys = new HashSet<>();
			for(final VepPrediction pred: this.vepPredictionParser.getPredictions(ctx)) {
				if(!accept(pred,ctx)) {
					if(isDebuggingVariant(ctx)) {
						LOG.info("VEP predictions  Ignoring "+shortName(ctx)+".");
					}
					continue;
				} else
				{
					if(isDebuggingVariant(ctx)) {
						LOG.info("VEP predictions  accepted "+shortName(ctx)+".");
					}
				}
				
				
				//ALL_NM && ALL_REFSEQ && ALL_ENST && ALL_TRANSCRIPTS
					{
					final String geneName =  pred.getSymbol();
					final String transcriptName =  pred.getFeature();
					
					
					if(!isEmpty(geneName) && !isEmpty(transcriptName)) {
						if(VcfBurdenSplitter2.this.enableAllNM && transcriptName.startsWith("NM_")) {
							keys.add(String.format("ALL_NM_%s_%s",ctx.getContig(),geneName));
							}
						if(VcfBurdenSplitter2.this.enableAllRefSeq && !transcriptName.startsWith("ENST"))
							{
							keys.add(String.format("ALL_REFSEQ_%s_%s",ctx.getContig(),geneName));
							}
						if(VcfBurdenSplitter2.this.enableAllEnst && transcriptName.startsWith("ENST"))
							{
							keys.add(String.format("ALL_ENST_%s_%s",ctx.getContig(),geneName));
							}
						}
					
					
					if(!isEmpty(geneName) && VcfBurdenSplitter2.this.enableAllGenes)
						{
						keys.add(String.format("ALL_GENES_%s_%s",ctx.getContig(),geneName));
						}
					if(!isEmpty(transcriptName) && VcfBurdenSplitter2.this.enableAllTranscript)
						{
						String k= String.format("ALL_TRANSCRIPTS_%s_%s",ctx.getContig(),transcriptName);
						if( !isEmpty(geneName) ) k+="_"+geneName;
						keys.add(k);
						}

					
					}
				
				String s;
				if(VcfBurdenSplitter2.this.enableVepHgnc) {
					s= pred.getHGNC();
					if(!isEmpty(s)) {
						keys.add(String.format("HGNC_%s_%s",ctx.getContig(),s));
						}
					}
				
				if(VcfBurdenSplitter2.this.enableVepEnsg) {
					s= pred.getEnsemblGene();
					if(!isEmpty(s)) {
						keys.add(String.format("ENSG_%s_%s",ctx.getContig(),s));
						}
					}
				/* same as feature 
				s= pred.getEnsemblTranscript();
				if(!isEmpty(s)) {
					keys.add(String.format("ENST_%s_%s",ctx.getContig(),s));
					}*/
				
				if(VcfBurdenSplitter2.this.enableVepFeature) {
					s= pred.getFeature();
					if(!isEmpty(s)) {
						keys.add(String.format("FEATURE_%s_%s",ctx.getContig(),s));
						
						if(VcfBurdenSplitter2.this.enableVepRefSeq && (s.startsWith("XM_") || s.startsWith("NM_")))
							{
							keys.add(String.format("REFSEQ_%s_%s",ctx.getContig(),s));
							}
						else if(VcfBurdenSplitter2.this.enableVepEnst && s.startsWith("ENST_"))
							{
							keys.add(String.format("ENST_%s_%s",ctx.getContig(),s));
							}
						}
					}
				
				if(VcfBurdenSplitter2.this.enableVepSymbol) {
					s= pred.getSymbol();
					if(!isEmpty(s)) {
						keys.add(String.format("SYMBOL_%s_%s",ctx.getContig(),s));
						}
					}
				
				if(VcfBurdenSplitter2.this.enableVepEnsp) {
					s= pred.getENSP();
					if(!isEmpty(s)) {
						keys.add(String.format("ENSP_%s_%s",ctx.getContig(),s));
						}
					}
				}
			return keys;
			}
		}
	
	private abstract class AbstractSoVepSplitter extends VepSplitter {
		final Set<SequenceOntologyTree.Term> acns;
		AbstractSoVepSplitter(final String acn_list[])
			{
			final SequenceOntologyTree soTree = SequenceOntologyTree.getInstance();
			this.acns = new HashSet<>(acn_list.length);
			for(final String soacn:acn_list)
				{
				final SequenceOntologyTree.Term  tacn = soTree.getTermByAcn(soacn);
				if(tacn==null)
					{
					throw new NullPointerException("tacn == null pour "+acns);
					}
				acns.addAll(tacn.getAllDescendants());
				}
			}
		@Override
		public boolean accept(final VepPrediction pred,final VariantContext origin) {
			for(final SequenceOntologyTree.Term so:pred.getSOTerms())
				{
				if(acns.contains(so))
					{
					if(isDebuggingVariant(origin)) {
						LOG.info("accepting variant "+shortName(origin)+" because SO-TERM "+so+" is in "+this.acns);
						}
					return true;
					}
				}
			if(isDebuggingVariant(origin)) {
				LOG.info("I don't accept variant "+shortName(origin)+" "+pred+" because SO-TERM "+pred.getSOTerms()+" is not in "+this.acns);
				}
			return false;
			}
		}
	
	private class SoVepSplitter extends AbstractSoVepSplitter {
		SoVepSplitter() {
			super(new String[]{
					"SO:0001893",  "SO:0001574",  "SO:0001575", 
					"SO:0001587",  "SO:0001589",  "SO:0001578", 
					"SO:0002012",  "SO:0001889",  "SO:0001821", 
					"SO:0001822",  "SO:0001583",  "SO:0001818"
					});
			}
		@Override
		public String getName() {
			return "vepso";
			}
		@Override
		public String getDescription() {
			return "Vep Sequence Ontology http://sequenceontology.org/ . Terms: "+ this.acns;
			}
	}
	
	private class HighDamageVepSplitter extends AbstractSoVepSplitter {
		HighDamageVepSplitter() {
			super(new String[]{
				"SO:0001893",  "SO:0001574",  "SO:0001575",
				"SO:0001587",  "SO:0001589",  "SO:0001578", 
				"SO:0002012",  "SO:0001889"
				});
		}
		@Override
		public String getName() {
			return "vephd";
			}
				
		@Override
		public String getDescription() {
			return "Vep Sequence Ontology http://sequenceontology.org/ . Terms:  High Damage. Terms: "+ this.acns;
			}
	}

	
	/** added for Matilde Nov 25, 2016, because all variants are FILTERed */
	private class ZeroVepFilter  extends AbstractSoVepSplitter {
		ZeroVepFilter() {
			super(new String[0]);
			}
		@Override
		public String getName() {
			return "vep0";
			}
		@Override
		public boolean accept(final VepPrediction pred,VariantContext ctx) {
			return true;
		}
		@Override
		public String getDescription() {
			return "Any Sequence Ontology term";
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
