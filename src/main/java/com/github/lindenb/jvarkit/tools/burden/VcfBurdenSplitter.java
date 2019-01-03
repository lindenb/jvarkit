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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC

### Description

This tools reads a VCF  (it should be sorted on chrom/POS and annotated with Ensembl Variation Predictor) and split data into genomic area of interest (gene, transcripts...).
For each area, a small VCF is produced and a Fished test is computed.
The final output is a set of concatenated VCF files. You could insert in a database using VcfDerby01


END_DOC
*/

@Program(name="vcfburdensplitter",description="Split VCF Using a Burden Splitter (by gene, transcript, etc..)")
public class VcfBurdenSplitter
	extends Launcher
	{
	
	private static final Logger LOG = Logger.build(VcfBurdenSplitter.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-ls","--listsplitters"},description="List available splitters and exit",help=true)
	private boolean listSplitter = false;

	@Parameter(names={"-sp","--splitterName"},description="Splitter Name")
	private String splitterName = "vepso";

	@Parameter(names={"-if","--ignorefilter"},description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
	private boolean acceptFiltered = false;

	@Parameter(names={"-gh","--galaxyhtml"},description="When used with galaxy, the files will be expanded in that path.")
	private String galaxyHtmlPath = "";

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
	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection=new WritingSortingCollection();
	
	//public for knime
	public static final String VCF_HEADER_SPLITKEY="VCFBurdenSplitName";
	
	
	/** getter for enableVepFeature */
	private boolean isEnableVepFeature()
		{
		return this.enableVepFeature;
		}

	/** getter for enableAllNM */
	private boolean isEnableAllNM()
		{
		return this.enableAllNM;
		}

	/** getter for enableAllGenes */
	private boolean isEnableAllGenes()
		{
		return this.enableAllGenes;
		}
	
	/** getter for enableVepEnsg */
	private boolean isEnableVepEnsg()
		{
		return this.enableVepEnsg;
		}
	
	/** getter for enableAllRefSeq */
	private boolean isEnableAllRefSeq()
		{
		return this.enableAllRefSeq;
		}
	
	/** getter for enableAllEnst */
	private boolean isEnableAllEnst()
		{
		return this.enableAllEnst;
		}

	/** getter for enableAllTranscript */
	private boolean isEnableAllTranscript()
		{
		return this.enableAllTranscript;
		}
	
	/** getter for enableVepHgnc */
	private boolean isEnableVepHgnc()
		{
		return this.enableVepHgnc;
		}

	/** getter for enableVepRefSeq */
	private boolean isEnableVepRefSeq()
		{
		return this.enableVepRefSeq;
		}
	
	
	/** getter for enableVepEnst */
	private boolean isEnableVepEnst()
		{
		return this.enableVepEnst;
		}
	
	
	/** getter for enableVepSymbol */
	private boolean isEnableVepSymbol()
		{
		return this.enableVepSymbol;
		}
	
	/** getter for enableVepEnsp */
	private boolean isEnableVepEnsp()
		{
		return this.enableVepEnsp;
		}

	
	private static String shortName(final VariantContext ctx) {
		return ctx.getContig()+":"+ctx.getStart()+":"+ctx.getID()+":"+ctx.getAlleles();
	}
	private boolean isDebuggingVariant(final VariantContext ctx) {
		return false;
	}
	
	private static class KeyAndLine {
		final String key;
		final String ctx;
		
		//volatile members used for fast comparaison
		String contig=null;
		int pos=-1;
		Allele ref=null;
		Allele alt=null;
		
		KeyAndLine(final String key,final String ctx) {
			this.key = key;
			this.ctx = ctx;
		}
		
		void build(final Pattern tab) {
			if(contig!=null) return;
			
			final String tokens1[] =  tab.split(this.ctx,6);
			this.contig=tokens1[0];
			this.pos=Integer.parseInt(tokens1[1]);
			this.ref = Allele.create(tokens1[3],true);
			this.alt = Allele.create(tokens1[4],false);
		}
		
	}
	
	private static class KeyAndLineComparator
		implements Comparator<KeyAndLine>
		{
		final Pattern tab = Pattern.compile("[\t]");
		@Override
		public int compare(final KeyAndLine o1, final KeyAndLine o2) {
			int i = o1.key.compareTo(o2.key);
			if(i!=0) return i;
			o1.build(tab);
			o2.build(tab);
			
			if(!o1.contig.equals(o2.contig)) {
				throw new IllegalStateException("not same contig???");
			}
			i =o1.pos - o2.pos;
			if(i!=0) return i;
			i = o1.ref.compareTo(o2.ref);
			if(i!=0) return i;
			return o1.alt.compareTo(o2.alt);
			}
		}
	
	private static class KeyAndLineCodec extends AbstractDataCodec<KeyAndLine>
		{
		@Override
		public KeyAndLine decode(final DataInputStream dis) throws IOException {
			String k;
			try {
				k=dis.readUTF();
			} catch(IOException err) { return null;}
			final String v = AbstractDataCodec.readString(dis);
			return new KeyAndLine(k, v);
		}
		@Override
		public void encode(final DataOutputStream dos, final KeyAndLine object) throws IOException {
			dos.writeUTF(object.key);
			AbstractDataCodec.writeString(dos, object.ctx);
			}
		@Override
		public AbstractDataCodec<KeyAndLine> clone() {
			return new KeyAndLineCodec();
			}
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
						if(isEnableAllNM() && transcriptName.startsWith("NM_")) {
							keys.add(String.format("ALL_NM_%s_%s",ctx.getContig(),geneName));
							}
						if(isEnableAllRefSeq() && !transcriptName.startsWith("ENST"))
							{
							keys.add(String.format("ALL_REFSEQ_%s_%s",ctx.getContig(),geneName));
							}
						if(isEnableAllEnst() && transcriptName.startsWith("ENST"))
							{
							keys.add(String.format("ALL_ENST_%s_%s",ctx.getContig(),geneName));
							}
						}
					
					
					if(!isEmpty(geneName) && isEnableAllGenes())
						{
						keys.add(String.format("ALL_GENES_%s_%s",ctx.getContig(),geneName));
						}
					if(!isEmpty(transcriptName) && isEnableAllTranscript())
						{
						String k= String.format("ALL_TRANSCRIPTS_%s_%s",ctx.getContig(),transcriptName);
						if( !isEmpty(geneName) ) k+="_"+geneName;
						keys.add(k);
						}

					
					}
				
				String s;
				if(isEnableVepHgnc()) {
					s= pred.getHGNC();
					if(!isEmpty(s)) {
						keys.add(String.format("HGNC_%s_%s",ctx.getContig(),s));
						}
					}
				
				if(isEnableVepEnsg()) {
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
				
				if(isEnableVepFeature()) {
					s= pred.getFeature();
					if(!isEmpty(s)) {
						keys.add(String.format("FEATURE_%s_%s",ctx.getContig(),s));
						
						if(isEnableVepRefSeq() && (s.startsWith("XM_") || s.startsWith("NM_")))
							{
							keys.add(String.format("REFSEQ_%s_%s",ctx.getContig(),s));
							}
						else if(isEnableVepEnst() && s.startsWith("ENST_"))
							{
							keys.add(String.format("ENST_%s_%s",ctx.getContig(),s));
							}
						}
					}
				
				if(isEnableVepSymbol()) {
					s= pred.getSymbol();
					if(!isEmpty(s)) {
						keys.add(String.format("SYMBOL_%s_%s",ctx.getContig(),s));
						}
					}
				
				if(isEnableVepEnsp()) {
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
	
	public VcfBurdenSplitter()
		{
		}
	
	
	
	@Override
	protected int doVcfToVcf(String inputName, File outorNull) {
		SortingCollection<KeyAndLine> sortingcollection=null;
		BufferedReader in = null;
		CloseableIterator<KeyAndLine> iter=null;
		PrintStream pw = null;
		PrintWriter allDiscardedLog = null;
		try {
			in = inputName==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(inputName)
					;
			
			if( this.allFilteredFileOut!=null) {
				allDiscardedLog = IOUtils.openFileForPrintWriter(this.allFilteredFileOut);
			}
					
			final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(in);
				
			/** find splitter by name */
			Splitter splitter = null;
			for(final Splitter s: this.splitters)
				{
				if(this.splitterName.equals(s.getName())) {
					splitter=s;
					break;
					}
				}
			if(splitter==null) {
				return wrapException("Cannot find a splitter named "+this.splitterName);
			}
			splitter.initialize(cah.header);
			LOG.info("splitter is "+splitter);
			
			
			pw= super.openFileOrStdoutAsPrintStream(outorNull);
			
			// read variants
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(cah.header);
			String prev_contig=null;
			for(;;)
				{
				final String line=in.readLine();
				final VariantContext variant = (line==null?null:progess.watch(cah.codec.decode(line)));
						
				if(variant==null || !variant.getContig().equals(prev_contig)) {
					if(sortingcollection!=null) {
						sortingcollection.doneAdding();
						iter = sortingcollection.iterator();
						LOG.info("dumping data for CONTIG: \""+prev_contig+"\"");
						
						final EqualRangeIterator<KeyAndLine> eqiter = new EqualRangeIterator<>(iter, new Comparator<KeyAndLine>() {
							@Override
							public int compare(final KeyAndLine o1, final KeyAndLine o2) {
								return o1.key.compareTo(o2.key);
								}
							});
						while(eqiter.hasNext())
							{
							final List<KeyAndLine> buffer = eqiter.next();
							
							final KeyAndLine first =  buffer.get(0);
							LOG.info(first.key);
							
							final List<VariantContext> variants = new ArrayList<>(buffer.size());
							boolean has_only_filtered=true;
							for(final KeyAndLine kal:buffer) {
								final VariantContext ctx = cah.codec.decode(kal.ctx);
								variants.add(ctx);
								
								if(isDebuggingVariant(ctx)) {
									LOG.info("Adding variant to list for key "+kal.key+" "+shortName(ctx));
									}
								
								if(!ctx.getContig().equals(prev_contig)) {
									eqiter.close();
									return wrapException("illegal state");
									}
								if(!ctx.isFiltered() || this.acceptFiltered) {
									has_only_filtered=false;
									//break; NOOOONNN !!!
									}
								}
							
							// all ctx are filtered			
							if(has_only_filtered)  {
								LOG.warn("ALL IS FILTERED in "+first.key);
								if( allDiscardedLog!=null) {
									for(final VariantContext ctx:variants) {
										if(isDebuggingVariant(ctx)) {
											LOG.info("Variant "+shortName(ctx)+" is part of never filtered for "+first.key);
											}
										
										allDiscardedLog.println(String.join("\t",
												first.key,
												ctx.getContig(),
												String.valueOf(ctx.getStart()),
												ctx.getReference().getDisplayString(),
												ctx.getAlternateAllele(0).getDisplayString(),
												String.valueOf(ctx.getFilters())
												));
										}
									}
								continue;
							}
							
							// save vcf file
							final VariantContextWriter out = VCFUtils.createVariantContextWriterToOutputStream(IOUtils.uncloseableOutputStream(pw));
							final VCFHeader header2=addMetaData(new VCFHeader(cah.header));
							header2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_SPLITKEY,first.key));
							
							out.writeHeader(header2);
							for(final VariantContext ctx:variants) {
								if(isDebuggingVariant(ctx))
									{
									LOG.info("saving variant "+shortName(ctx)+" to final output with key="+first.key);
									}
								out.add(ctx);
							}
							out.close();//yes because wrapped into IOUtils.encloseableOutputSream
							pw.flush();
							}
						eqiter.close();
						iter.close();iter=null;
						//dispose sorting collection
						sortingcollection.cleanup();
						sortingcollection=null;
						}
					//EOF met
					if(variant==null) break;
					prev_contig = variant.getContig();
					}
				
				if(sortingcollection==null) {
					/* create sorting collection for new contig */
					sortingcollection = SortingCollection.newInstance(
							KeyAndLine.class,
							new KeyAndLineCodec(),
							new KeyAndLineComparator(),
							this.writingSortingCollection.maxRecordsInRam,
							this.writingSortingCollection.getTmpPaths()
							);
					sortingcollection.setDestructiveIteration(true);
					}
				
				if( variant.getAlternateAlleles().size()!=1) {
					return wrapException("Expected only one allele per variant. Please use VcfMultiToOneAllele https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele.");
					}
								

				
				//no check for ctx.ifFiltered here, we do this later.
				for(final String key: splitter.keys(variant)) {
					if(isDebuggingVariant(variant)) {
						LOG.info("Adding variant with key "+key+" "+shortName(variant));
						}
					sortingcollection.add(new KeyAndLine(key, line));
					}
				}
			progess.finish();
			
			
			pw.flush();
			pw.close();pw=null;
			
			if(allDiscardedLog!=null)
				{
				allDiscardedLog.flush();
				allDiscardedLog.close();
				allDiscardedLog=null;
				}
			
			return RETURN_OK;
			}
		catch(final Exception err) 
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(iter);
			if(sortingcollection!=null) sortingcollection.cleanup();
			CloserUtil.close(in);
			CloserUtil.close(pw);
			CloserUtil.close(allDiscardedLog);
			}
		}
		@Override
		public int doWork(List<String> args) {
		
		if(!this.listSplitter) {
			if(this.outputFile==null) {
				LOG.warn("output file option -o was not be declared. Zip will be printed to stdout.");
			}
			else if(!outputFile.getName().endsWith(".zip")) {
				LOG.error("output file option -o is not be declared a en with .zip");
				return -1;
			}
			
			
			if(this.splitterName.isEmpty()) {
				return wrapException("splitter name undefined");
			}
			Splitter splitter=null;
			for(final Splitter s: this.splitters)
				{
				if(this.splitterName.equals(s.getName())) {
					splitter=s;
					break;
					}
				}
			if(splitter==null) {
			return wrapException("Cannot find a splitter named "+this.splitterName);
			}
		}
		if(this.listSplitter) {
			for(final Splitter splitter:this.splitters) {
				stdout().println(splitter.getName()+"\t"+splitter.getDescription());
	
			}
			return RETURN_OK;
			}
		
		
		return doVcfToVcf(args,this.outputFile);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenSplitter().instanceMainWithExit(args);
		}
	}
