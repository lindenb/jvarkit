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
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
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
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.DelegateVariantContextWriter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.AnnPredictionParser.AnnPrediction;

/**

 * @author lindenb
 * 
 BEGIN_DOC
 
 
 
 END_DOC
 
 */
@Program(name="vcfburdensplitter2",description="new version",keywords={"vcf","burden","gene","vep","snpeff","prediction"})
public class VcfBurdenSplitter2 
	{
	//public for knime
	public static final String DEFAULT_VCF_HEADER_SPLITKEY="VCFBurdenSplitName";

	private static final Logger LOG = Logger.build().
			prefix("splitter2").
			make();

	@Parameter(names={"-if","--ignorefilter"},description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
	private boolean acceptFiltered = false;

	@Parameter(names="--maxRecordsInRam",description="Max records in RAM")
	private int maxRecordsInRam=50000;
			
	@Parameter(names="--tmpDir",description="Temporary directory")
	private File tmpDir = IOUtils.getDefaultTmpDir();
	
	@Parameter(names={"-m","--manifestFile"},description="Manifest File")
	private File manifestFile = null;

	@Parameter(names={"-t","--tag"},description="Split Key")
	private String splitInfoKey = DEFAULT_VCF_HEADER_SPLITKEY;

	
	
	public VariantContextWriter open(final VariantContextWriter delegate) {
		final MyWriter w = new MyWriter(delegate);
		w.manifestFile=this.manifestFile;
		w.acceptFiltered=this.acceptFiltered;
		w.tmpDir=this.tmpDir;
		w.maxRecordsInRam=this.maxRecordsInRam;
		w.splitInfoKey=this.splitInfoKey; 
		return w;
		}
	
	
	private  class MyWriter extends DelegateVariantContextWriter
		{
		private String prev_contig=null;
		private SortingCollection<Interval> sortingcollection=null;
		private PrintWriter manifestWriter=null;
		private File manifestFile = null;
		private boolean acceptFiltered = false;
		private File tmpDir=null;
		private int maxRecordsInRam;
		private String splitInfoKey;
		private AnnPredictionParser annPredictionParser = null;
		private VepPredictionParser vepPredictionParser = null;

		MyWriter(final VariantContextWriter w)
			{
			super(w);
			}
		
		@Override
		public void writeHeader(final VCFHeader header) {
			final VCFHeader header2= new VCFHeader(header);
			
			this.annPredictionParser = new AnnPredictionParserFactory(header).get();
			this.vepPredictionParser = new VepPredictionParserFactory(header).get();
			
			this.prev_contig=null;
			if(this.manifestFile==null) {
				LOG.warning("Manifest file is undefined");
				this.manifestWriter=new PrintWriter(new NullOuputStream());
				}
			else
				{
				try {
					this.manifestWriter=new PrintWriter(this.manifestFile);
				} catch (final FileNotFoundException e) {
					throw new RuntimeIOException(e);
				}
				}
			
			
			
			header2.addMetaDataLine(
				new VCFInfoHeaderLine(
						this.splitInfoKey,
						VCFHeaderLineCount.UNBOUNDED,
						VCFHeaderLineType.String,
						"Split Names"
						));
			super.writeHeader(header2);
			}
		
		@Override
		public void add(final VariantContext ctx) {
			if(prev_contig==null || !ctx.getContig().equals(prev_contig)) {
				dump();
				prev_contig=ctx.getContig();
				}
			if(ctx.isFiltered() && !this.acceptFiltered)
				{
				//add to delegate without registering the SPLIT name
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
						this.maxRecordsInRam,
						this.tmpDir.toPath()
						);
				this.sortingcollection.setDestructiveIteration(true);
				}
			for(final String spltiName:splitNames)
				{
				this.sortingcollection.add(
						new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd(),false,spltiName));
				}
			
			
			super.add(new VariantContextBuilder(ctx).attribute(
					this.splitInfoKey,
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
				this.manifestWriter=null;
				}
			super.close();
			}
		
		
		private Set<String> getSplitNamesFor(final VariantContext ctx){
			final Set<String> keys = new HashSet<>();
			for(final VepPrediction pred: this.vepPredictionParser.getPredictions(ctx)) {
				keys.addAll(pred.getGeneKeys().stream().map(S->ctx.getContig()+"_"+S).collect(Collectors.toSet()));
				}
			for(final AnnPrediction pred: this.annPredictionParser.getPredictions(ctx)) {
				keys.addAll(pred.getGeneKeys().stream().map(S->ctx.getContig()+"_"+S).collect(Collectors.toSet()));
				}	
			/* replace . by _ so we don't have problems with regex later */
			return keys.stream().
					map(S->S.replace('.', '_').replace('-', '_')).
					collect(Collectors.toSet());
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
			if(this.manifestWriter.checkError()) {
				LOG.warn("I/O error when writing manifest");
			}
			eqiter.close();
			iter.close();iter=null;
			//dispose sorting collection
			sortingcollection.cleanup();
			sortingcollection=null;
			}

		
		}
	
	

	
	
	
	
	
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

	protected boolean isDebuggingVariant(final VariantContext ctx) {
		return false;
	}
	
	protected String shortName(final VariantContext ctx) {
	return ctx.getContig()+":"+ctx.getStart()+":"+ctx.getAlleles();	
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
		@Parameter(names={"-o","--out"},description="Vcf output.")
		private VariantContextWriter output=new com.github.lindenb.jvarkit.util.jcommander.Launcher.VcfWriterOnDemand();
		
		Launcher() {
			
		
			}
		@Override
		public int doWork(final List<String> args) {
			VariantContextWriter w=null;
			VCFIterator in=null;
			try
				{
				in = VCFUtils.createVCFIterator(super.oneFileOrNull(args));
				w= this.instance.open(output);
				VCFUtils.copyHeaderAndVariantsTo(in, w);
				w.close();
				return 0;
				}
			catch(Exception err)
				{
				LOG.fatal(err);
				return -1;
				}
			finally
				{
				CloserUtil.close(w);
				CloserUtil.close(in);
				}
			}
		}
	
	
	public static void main(String[] args)
		{
		new Launcher().instanceMainWithExit(args);
		}
	}
