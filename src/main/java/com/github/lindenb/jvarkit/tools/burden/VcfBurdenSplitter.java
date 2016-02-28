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
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.util.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;


public class VcfBurdenSplitter
	extends AbstractVcfBurdenSplitter
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenSplitter.class);
	
	private enum SuperVariant
		{
		SV0,AT_LEAST_ONE_VARIANT
		}
	
	private static class KeyAndLine {
		final String key;
		final String ctx;
		KeyAndLine(String key,String ctx) {
			this.key = key;
			this.ctx = ctx;
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
			final String tokens1[] = this.tab.split(o1.ctx,6);
			final String tokens2[] = this.tab.split(o2.ctx,6);
			if(!tokens1[0].equals(tokens2[0])) {
				throw new IllegalStateException("not same contig???");
			}
			i = Integer.parseInt(tokens1[1]) - Integer.parseInt(tokens2[1]);
			if(i!=0) return i;
			i = Allele.create(tokens1[3],true).compareTo( Allele.create(tokens2[3],true));
			if(i!=0) return i;
			return Allele.create(tokens1[4],false).compareTo( Allele.create(tokens2[4],false));
			}
		}
	private static class KeyAndLineCodec extends AbstractDataCodec<KeyAndLine>
		{
		@Override
		public KeyAndLine decode(DataInputStream dis) throws IOException {
			String k;
			try {
				k=dis.readUTF();
			} catch(IOException err) { return null;}
			final String v = AbstractDataCodec.readString(dis);
			return new KeyAndLine(k, v);
		}
		@Override
		public void encode(DataOutputStream dos, KeyAndLine object) throws IOException {
			dos.writeUTF(object.key);
			AbstractDataCodec.writeString(dos, object.ctx);
			}
		@Override
		public AbstractDataCodec<KeyAndLine> clone() {
			return new KeyAndLineCodec();
			}
		}

	
	private abstract class Splitter {
		public abstract Set<String> keys(final VariantContext ctx);
		public abstract String getName();
		public abstract String getDescription();
		public void initialize(final VCFHeader header) {}
		}
	
		private class VepSplitter extends Splitter {
		private VepPredictionParser vepPredictionParser=null;
		VepSplitter() {}
		
		
		
		@Override public String getName() { return "vep";}
		@Override public String getDescription() { return "Ensembl Variant Effect Prediction";}
		public void initialize(final VCFHeader header) {
			this.vepPredictionParser = new VepPredictionParser(header);
			}
		public boolean accept(VepPrediction pred) {
			return true;
		}
		private boolean isEmpty(final String s) {
			return s==null || s.trim().isEmpty();
		}
		@Override
		public Set<String> keys(final VariantContext ctx) {
			final Set<String> keys = new HashSet<>();
			for(final VepPrediction pred: this.vepPredictionParser.getPredictions(ctx)) {
				if(!accept(pred)) continue;
				String s= pred.getHGNC();
				if(!isEmpty(s)) {
					keys.add(String.format("HGNC_%s_%s",ctx.getContig(),s));
					}
				s= pred.getEnsemblGene();
				if(!isEmpty(s)) {
					keys.add(String.format("ENSG_%s_%s",ctx.getContig(),s));
					}
				s= pred.getEnsemblTranscript();
				if(!isEmpty(s)) {
					keys.add(String.format("ENST_%s_%s",ctx.getContig(),s));
					}
				s= pred.getFeature();
				if(!isEmpty(s)) {
					keys.add(String.format("FEATURE_%s_%s",ctx.getContig(),s));
					}
				s= pred.getSymbol();
				if(!isEmpty(s)) {
					keys.add(String.format("SYMBOL_%s_%s",ctx.getContig(),s));
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
		public boolean accept(final VepPrediction pred) {
			boolean ok=false;
			for(SequenceOntologyTree.Term so:pred.getSOTerms())
				{
				if(acns.contains(so))
					{
					ok=true;
					}
				}
			return ok;
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
			return "Vep SO. Terms: "+ this.acns;
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
			return "Vep High Damage. Terms: "+ this.acns;
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
	
	private final Splitter splitters[]= new Splitter[]{
			new VepSplitter(),
			new SoVepSplitter(),
			new HighDamageVepSplitter(),
			new SlidingWindowSplitter(1000, 500),
			new SlidingWindowSplitter(1000, 300),
		};
	
	public VcfBurdenSplitter()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.casesFile==null || !super.casesFile.exists()) {
		return wrapException("Undefined Case file option -"+OPTION_CASESFILE);
		}
		if(super.controlsFile==null || !super.controlsFile.exists()) {
		return wrapException("Undefined Control file option -"+OPTION_CONTROLSFILE);
		}
		if(super.outputFile==null || !outputFile.getName().endsWith(".zip")) {
			return wrapException("output file option -"+OPTION_OUTPUTFILE+" must be declared and must en with .zip");
		}
		
		if(!super.listSplitter) {
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
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			String inputName,
			VcfIterator in,
			VariantContextWriter out
			) throws IOException {
		throw new IllegalArgumentException("should be never called");
		}
	
	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName) throws Exception {
		SortingCollection<KeyAndLine> sortingcollection=null;
		BufferedReader in = null;
		FileOutputStream fos = null;
		ZipOutputStream zout=null;
		CloseableIterator<KeyAndLine> iter=null;
		PrintWriter pw = null;
		try {
			in = inputName==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(inputName)
					;
			
			final VCFUtils.CodecAndHeader cah = VCFUtils.parseHeader(in);
			final Set<Pedigree.Person> caseSamples = Pedigree.readPedigree(super.casesFile).getPersons();
			final Set<Pedigree.Person> controlSamples = Pedigree.readPedigree(super.controlsFile).getPersons();
			for(int pop=0;pop<2;++pop)
				{
				Iterator<Pedigree.Person> it = (pop==0?caseSamples:controlSamples).iterator();
				while(it.hasNext()) {
					Pedigree.Person p = it.next();
					if(!cah.header.getSampleNamesInOrder().contains(p.getId())) {
						LOG.warn("REMOVING "+p+" as it is not in vcf header");
						it.remove();
					}
				}
				}
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
			
			sortingcollection = SortingCollection.newInstance(
					KeyAndLine.class,
					new KeyAndLineCodec(),
					new KeyAndLineComparator(),
					super.maxRecordsInRam,
					super.getTmpDirectories()
					);
			sortingcollection.setDestructiveIteration(true);
			
			SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(cah.header.getSequenceDictionary());
			
			String line;
			while((line=in.readLine())!=null)
				{
				final VariantContext ctx = progess.watch(cah.codec.decode(line));
				if(	ctx.getAlternateAlleles().size()!=1) {
					return wrapException("Expected only one allele per variant. Please use ManyAlleletoOne.");
					}
				
				//no check for ctx.ifFiltered here, we do this later.
				
				
				for(final String key: splitter.keys(ctx)) {
					sortingcollection.add(new KeyAndLine(key, line));
					}
				}
			progess.finish();
			sortingcollection.doneAdding();
			
			LOG.info("creating zip "+super.outputFile);
			fos = new FileOutputStream(super.outputFile);
			zout = new ZipOutputStream(fos);
			

			final File tmpReportFile = File.createTempFile("_tmp.", ".txt", super.getTmpdir());
			tmpReportFile.deleteOnExit();
			pw = IOUtils.openFileForPrintWriter(tmpReportFile);
			pw.println("#chrom\tstart\tend\tkey\tFisher\tCount_Variants\tCount_non_filtered_Variants");
			
			iter = sortingcollection.iterator();
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
				String contig=null;
				int chromStart=Integer.MAX_VALUE;
				int chromEnd=0;
				int count_non_filtered=0;
				for(final KeyAndLine kal:buffer) {
					final VariantContext ctx = cah.codec.decode(kal.ctx);
					variants.add(ctx);
					if(!ctx.isFiltered() || super.acceptFiltered) {
						count_non_filtered++;
						}
					contig = ctx.getContig();
					chromStart = Math.min( chromStart , ctx.getStart() );
					chromEnd = Math.max( chromEnd , ctx.getEnd() );
					}
				
			// all ctx are filtered			
			if( count_non_filtered == 0)  continue;
			
			/* initialize supervariants */
			int count_case_sv0=0;
			int count_ctrl_sv0=0;
			int count_case_sv1=0;
			int count_ctrl_sv1=0;


			//loop over case control
			for(int pop=0;pop<2;++pop) {
				for(final Pedigree.Person person : (pop==0?caseSamples:controlSamples)) {
					SuperVariant superVariant = SuperVariant.SV0;
					
					for(final VariantContext ctx : variants)
						{
						if(ctx.isFiltered() && !super.acceptFiltered) continue;

						final Genotype g = ctx.getGenotype(person.getId());	
						if(g==null) continue;//not in vcf header
						final Allele alt = ctx.getAlternateAlleles().get(0);

						for(final Allele a:g.getAlleles())
							{
							if(a.equals(alt)) {
								superVariant = SuperVariant.AT_LEAST_ONE_VARIANT;
								break;
								}
							}//end of allele
						if(superVariant!=SuperVariant.SV0 ) break;
						}//end of variant
					if(superVariant==SuperVariant.SV0 ) {
						if(pop==0 ) count_case_sv0++;
						else count_ctrl_sv0++;
					} else // AT_LEAST_ONE_VARIANT 
						{
						if(pop==0 ) count_case_sv1++;
						else count_ctrl_sv1++;
						}
						
					}//end of case/ctrl
				}//end of pop
				
				
	
				final FisherExactTest fisher = FisherExactTest.compute(
						count_case_sv0, count_case_sv1,
						count_ctrl_sv0, count_ctrl_sv1
						);
				pw.println(
						contig+"\t"+
						(chromStart-1)+"\t"+//-1 for bed compatibility
						chromEnd+"\t"+
						first.key+"\t"+
						fisher.getAsDouble()+"\t"+
						variants.size()+"\t"+
						count_non_filtered
						);
				
				// save vcf file
				final ZipEntry ze = new ZipEntry(super.baseZipDir+"/"+first.key+".vcf");
				zout.putNextEntry(ze);
				final VariantContextWriter out = VCFUtils.createVariantContextWriterToOutputStream(zout);
				final VCFHeader header2=addMetaData(new VCFHeader(cah.header));
				header2.addMetaDataLine(new VCFHeaderLine("VCFBurdenFisher",
						String.valueOf(fisher.getAsDouble())));
				out.writeHeader(header2);
				for(final VariantContext ctx:variants) {
					out.add(ctx);
				}
				//out.close();//NON
				zout.closeEntry();
				}
			eqiter.close();
			iter.close();iter=null;
			
			progess.finish();
			
			LOG.info("saving report");
			pw.flush();
			pw.close();
			final ZipEntry entry = new ZipEntry(super.baseZipDir+"/fisher.tsv");
			zout.putNextEntry(entry);
			IOUtils.copyTo(tmpReportFile,zout);
			zout.closeEntry();
			
			zout.finish();
			zout.close();
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
			CloserUtil.close(fos);
			CloserUtil.close(pw);
			}
		}
		
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
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
