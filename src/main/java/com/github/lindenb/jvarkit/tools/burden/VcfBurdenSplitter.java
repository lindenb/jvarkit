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
import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


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
		}
	
	private class SlidingWindowSplitter extends Splitter {
		final int winsize;
		final int winshift;
		SlidingWindowSplitter(final VCFHeader header,int winsize,int winshift)
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
	}
	
	
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
			
			Splitter splitter = new SlidingWindowSplitter(cah.header, 1000, 500);
			
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
			pw.println("#key\tFisher\tSuperVariant");
			
			iter = sortingcollection.iterator();
			final List<KeyAndLine> buffer=new ArrayList<>();
			for(;;)
				{
				KeyAndLine curr = null;
				if( iter.hasNext() ) {
					curr = iter.next();
					}
				if(curr==null || (!buffer.isEmpty() && !buffer.get(0).key.equals(curr.key))) {
					if(!buffer.isEmpty()) {
						final KeyAndLine first =  buffer.get(0);
						LOG.info(first.key);
						final List<VariantContext> variants = new ArrayList<>(buffer.size());
						for(final KeyAndLine kal:buffer) {
							final VariantContext ctx = cah.codec.decode(kal.ctx);
							variants.add(ctx);
						}
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
						pw.println(first.key+"\t"+fisher.getAsDouble());
						
						
						ZipEntry ze = new ZipEntry(super.baseZipDir+"/"+first.key+".vcf");
						zout.putNextEntry(ze);
						VariantContextWriter out = VCFUtils.createVariantContextWriterToOutputStream(zout);
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
					if(curr==null) break;
					buffer.clear();
					}
				buffer.add(curr);
				}
			iter.close();iter=null;
			
			progess.finish();
			
			LOG.info("saving report");
			pw.flush();
			pw.close();
			ZipEntry entry = new ZipEntry(super.baseZipDir+"/fisher.tsv");
			zout.putNextEntry(entry);
			IOUtils.copyTo(tmpReportFile,zout);
			zout.closeEntry();
			
			zout.finish();
			zout.close();
			return RETURN_OK;
			}
		catch(Exception err) 
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
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenSplitter().instanceMainWithExit(args);
		}
	}
