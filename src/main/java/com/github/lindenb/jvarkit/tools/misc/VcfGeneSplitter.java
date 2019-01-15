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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser.VepPrediction;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**

BEGIN_DOC

### Example

```

$ gunzip -c input.vcf.gz |\
	java -jar dist/vcfgenesplitter.jar -tmpdir . -o out.zip -maxRecordsInRam 5000 -zipdir BASE

```

END_DOC
*/
@Program(name="vcfgenesplitter",description="Split VCF+VEP by gene.")
public class VcfGeneSplitter
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfGeneSplitter.class).make();
	
	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;


	@Parameter(names={"-zipdir","--zipdir"},description="Base zip entry")
	private String baseZipDir = "splitbygene";

	@Parameter(names={"-enst","--noENST"},description="Ignore ENST entries")
	private boolean ignoreEnst = false;

	@Parameter(names={"-hgns","--noHGNC"},description="Ignore HGNC entries")
	private boolean ignoreHgnc = false;

	@Parameter(names={"-feature","--noFeature"},description="Ignore FEATURE entries")
	private boolean ignoreFeature = false;

	@Parameter(names={"-ensg","--noENSG"},description="Ignore ENSG entries")
	private boolean ignoreEnsg = false;

	@Parameter(names={"-refseq","--noREFSEQ"},description="Ignore ENSG entries")
	private boolean ignoreRefSeq = false;

	@Parameter(names={"-symbol","--noSymbol"},description="Ignore SYMBOL entries")
	private boolean ignoreSymbol = false;

	@Parameter(names={"-ensp","--noENSP"},description="Ignore ENSP entries")
	private boolean ignoreENSP = false;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();
	
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
			final String tokens1[] = this.tab.split(o1.ctx,5);
			final String tokens2[] = this.tab.split(o2.ctx,5);
			if(!tokens1[0].equals(tokens2[0])) {
				throw new IllegalStateException("not same contig???");
			}
			i = Integer.parseInt(tokens1[1]) - Integer.parseInt(tokens2[1]);
			if(i!=0) return i;
			i = Allele.create(tokens1[3],true).compareTo( Allele.create(tokens2[3],true));
			return i;
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

	
	private static boolean isEmpty(final String s) {
		return s==null || s.trim().isEmpty();
	}

	private Set<String> getVariantKeys(final VepPredictionParser vepPredictionParser,final VariantContext ctx)
			{
			final Set<String> keys = new HashSet<>();
			for(final VepPrediction pred: vepPredictionParser.getPredictions(ctx)) {
				String s= pred.getHGNC();
				if(!isEmpty(s) && !this.ignoreHgnc) {
					keys.add(String.format("HGNC_%s_%s",ctx.getContig(),s));
					}
				s= pred.getEnsemblGene();
				if(!isEmpty(s) && !this.ignoreEnsg) {
					keys.add(String.format("ENSG_%s_%s",ctx.getContig(),s));
					}
				/* same as feature 
				s= pred.getEnsemblTranscript();
				if(!isEmpty(s)) {
					keys.add(String.format("ENST_%s_%s",ctx.getContig(),s));
					}*/
					
				s= pred.getFeature();
				if(!isEmpty(s)  && !this.ignoreFeature) {
					keys.add(String.format("FEATURE_%s_%s",ctx.getContig(),s));
					
					if((s.startsWith("XM_") || s.startsWith("NM_")) &&  !this.ignoreRefSeq)
						{
						keys.add(String.format("REFSEQ_%s_%s",ctx.getContig(),s));
						}
					else if(s.startsWith("ENST_")&&  !this.ignoreEnst)
						{
						keys.add(String.format("ENST_%s_%s",ctx.getContig(),s));
						}
					}

				s= pred.getSymbol();
				if(!isEmpty(s) && !this.ignoreSymbol) {
					keys.add(String.format("SYMBOL_%s_%s",ctx.getContig(),s));
					}
				s= pred.getENSP();
				if(!isEmpty(s) && !this.ignoreENSP) {
					keys.add(String.format("ENSP_%s_%s",ctx.getContig(),s));
					}
				}
			return keys;
			}
		
	
	public VcfGeneSplitter()
		{
		}
	
	
	
	
	@Override
	protected int doVcfToVcf(String inputName, File outputFile) {
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
			/** find splitter by name */
			
			final VepPredictionParser vepPredictionParser=new VepPredictionParserFactory().header(cah.header).get();
			sortingcollection = SortingCollection.newInstance(
					KeyAndLine.class,
					new KeyAndLineCodec(),
					new KeyAndLineComparator(),
					this.writingSortingCollection.getMaxRecordsInRam(),
					this.writingSortingCollection.getTmpPaths()
					);
			sortingcollection.setDestructiveIteration(true);
			
			// read variants
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(cah.header);
			String line;
			
			while((line=in.readLine())!=null)
				{
				final VariantContext ctx = progess.watch(cah.codec.decode(line));
				//no check for ctx.ifFiltered here, we do this later.
				for(final String key: this.getVariantKeys(vepPredictionParser,ctx)) {
					sortingcollection.add(new KeyAndLine(key, line));
					}
				}
			progess.finish();
			sortingcollection.doneAdding();
			
			LOG.info("creating zip "+ outputFile);
			fos = new FileOutputStream(outputFile);
			zout = new ZipOutputStream(fos);
			

			final File tmpReportFile = File.createTempFile("_tmp.", ".txt", writingSortingCollection.getTmpDirectories().get(0));
			tmpReportFile.deleteOnExit();
			pw = IOUtils.openFileForPrintWriter(tmpReportFile);
			pw.println("#chrom\tstart\tend\tkey\tCount_Variants");
			
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
				for(final KeyAndLine kal:buffer) {
					final VariantContext ctx = cah.codec.decode(kal.ctx);
					variants.add(ctx);
					contig = ctx.getContig();
					chromStart = Math.min( chromStart , ctx.getStart() );
					chromEnd = Math.max( chromEnd , ctx.getEnd() );
					}
				

				
				
				pw.println(
						contig+"\t"+
						(chromStart-1)+"\t"+//-1 for bed compatibility
						chromEnd+"\t"+
						first.key+"\t"+
						variants.size()
						);
				
				// save vcf file
				final ZipEntry ze = new ZipEntry(this.baseZipDir+"/"+first.key+".vcf");
				zout.putNextEntry(ze);
				final VariantContextWriter out = VCFUtils.createVariantContextWriterToOutputStream(IOUtils.uncloseableOutputStream(zout));
				final VCFHeader header2=addMetaData(new VCFHeader(cah.header));
				header2.addMetaDataLine(new VCFHeaderLine("VcfGeneSplitter.Name",
						String.valueOf(first.key)));
				out.writeHeader(header2);
				for(final VariantContext ctx:variants) {
					out.add(ctx);
				}
				out.close();//yes because wrapped into IOUtils.encloseableOutputSream
				zout.closeEntry();
				}
			eqiter.close();
			iter.close();iter=null;
			
			progess.finish();
			
			LOG.info("saving report");
			pw.flush();
			pw.close();
			final ZipEntry entry = new ZipEntry(this.baseZipDir+"/manifest.bed");
			zout.putNextEntry(entry);
			IOUtils.copyTo(tmpReportFile,zout);
			zout.closeEntry();
			
			zout.finish();
			zout.close();
			return RETURN_OK;
			}
		catch(final Exception err) 
			{
			LOG.error(err);
			return -1;
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
	public int doWork(List<String> args) {
		if(this.outputFile==null || !outputFile.getName().endsWith(".zip")) {
			LOG.error("output file option -o must be declared and must en with .zip");
			return -1;
			}	
		return doVcfToVcf(args, outputFile);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfGeneSplitter().instanceMainWithExit(args);
		}
	}
