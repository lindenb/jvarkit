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

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
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
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;


public class VcfBurdenSplitter
	extends AbstractVcfBurdenSplitter
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenSplitter.class);
	
	//public for knime
	public static final String VCF_HEADER_SPLITKEY="VCFBurdenSplitName";
	
	
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
			this.vepPredictionParser = new VepPredictionParser(header);
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
		public boolean accept(final VepPrediction pred) {
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
		SlidingWindowSplitter(int winsize,int winshift)
			{
			this.winsize = winsize;
			this.winshift =winshift;
			}
		@Override
		public Set<String> keys(final VariantContext ctx)
			{
			final Set<String> set= new HashSet<>();
			int x= (int)(ctx.getStart()/this.winsize);
			x*=winsize;
			while(x<=ctx.getStart() && ctx.getStart()<=x+this.winshift) {
				set.add(String.format("%s_%09d_%09d", ctx.getContig(),x,x+winshift));
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
		};
	
	public VcfBurdenSplitter()
		{
		}
	
	
	@Override
	public Collection<Throwable> initializeKnime() {
		
		if(!super.listSplitter) {
			if(super.outputFile==null) {
				LOG.warn("output file option -"+OPTION_OUTPUTFILE+" was not be declared. Zip will be printed to stdout.");
			}
			else if(!outputFile.getName().endsWith(".zip")) {
				return wrapException("output file option -"+OPTION_OUTPUTFILE+" is not be declared a en with .zip");
			}
			
			
			if(super.splitterName.isEmpty()) {
				return wrapException("splitter name undefined in option -"+OPTION_SPLITTERNAME);
			}
			Splitter splitter=null;
			for(final Splitter s: this.splitters)
				{
				if(super.splitterName.equals(s.getName())) {
					splitter=s;
					break;
					}
				}
			if(splitter==null) {
			return wrapException("Cannot find a splitter named "+super.splitterName);
			}
		}
		return super.initializeKnime();
	 	}
	
	
	//@Override
	public Collection<Throwable> doVcfToVcf(final String inputName) throws Exception {
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
			
			if( super.allFilteredFileOut!=null) {
				allDiscardedLog = IOUtils.openFileForPrintWriter(super.allFilteredFileOut);
			}
					
			final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(in);
				
			/** find splitter by name */
			Splitter splitter = null;
			for(final Splitter s: this.splitters)
				{
				if(super.splitterName.equals(s.getName())) {
					splitter=s;
					break;
					}
				}
			if(splitter==null) {
				return wrapException("Cannot find a splitter named "+super.splitterName);
			}
			splitter.initialize(cah.header);
			LOG.info("splitter is "+splitter);
			
			
			pw= super.openFileOrStdoutAsPrintStream();
			
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
								if(!ctx.isFiltered() || super.acceptFiltered) {
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
									LOG.info("saving variant "+shortName(ctx)+" to final output");
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
							super.maxRecordsInRam,
							super.getTmpDirectories()
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
	protected Collection<Throwable> call(final String inputName) throws Exception {
		if(super.listSplitter) {
			for(final Splitter splitter:this.splitters) {
				stdout().println(splitter.getName()+"\t"+splitter.getDescription());
	
			}
			return RETURN_OK;
		}
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenSplitter().instanceMainWithExit(args);
		}
	}
