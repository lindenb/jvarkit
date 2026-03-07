/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.rvtests;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.predictions.GeneExtractorFactory.GeneExtractor;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class RVTestsMakeSetFiles extends Launcher {
	private static final Logger LOG = Logger.of(RVTestsMakeSetFiles.class);

	@Parameter(names = "-o",description ="output directory",required = true)
	private Path outputDirectory = null;

	@ParametersDelegate
	private WritingSortingCollection writingSortingCollection = new WritingSortingCollection();


	
	private static class Variation {
		int tid;
		int start;
		int end;
		Allele ref;
		Allele alt;
		String key;

		
		int compareVariantOnKey(Variation v) {
			return this.key.compareTo(v.key);
			}
		int compareOnCoordinate(Variation v) {
			int i =Integer.compare(this.tid, v.tid);// CHROMOSOME
			if (i != 0) return i;
			i = Integer.compare(this.start, v.start);
			if (i != 0) return i;
			i = ref.compareTo(v.ref);
			if (i != 0) return i;
			i = alt.compareTo(v. alt);
			return i;
		 	}
		
		int compareVariantOnKeyAndPosition(Variation v) {
			int i = compareVariantOnKey(v);
			if (i != 0) return i;
			i =Integer.compare(this.tid, v.tid);// CHROMOSOME
			if (i != 0) return i;
			i = Integer.compare(this.start, v.start);// CHROMOSOME
			if (i != 0) return i;
			return Integer.compare(this.end, v.end);
		 	}
		
		@Override
		public String toString() {
			return tid+":"+start+":"+ref.getDisplayString()+":"+alt.getDisplayString();
			}
		}



	private static class RowCodec extends AbstractDataCodec<Variation> {
		@Override
		public void encode(DataOutputStream dos, final Variation variant) throws IOException {
			dos.writeInt(variant.tid);
			dos.writeInt(variant.start);
			dos.writeInt(variant.end);
			dos.writeUTF(variant.ref.getDisplayString());
			dos.writeUTF(variant.alt.getDisplayString());
			dos.writeUTF(variant.key);
			}
		
		@Override
		public Variation decode(DataInputStream dis) throws IOException {
			final Variation v = new Variation();
			try {
				v.tid = dis.readInt();
			} catch (EOFException err) {
				return null;
			}
			v.start = dis.readInt();
			v.end = dis.readInt();
			v.ref = Allele.create( dis.readUTF(),true);
			v.alt = Allele.create( dis.readUTF(),false);
			v.key = dis.readUTF();
			return v;
			}

		@Override
		public RowCodec clone() {
			return new RowCodec();
		}
	}
	
	protected String fixContig(final String ctg) {
		if(ctg.startsWith("chr")) return ctg.substring(3);
		return ctg;
		}

	
	private abstract class Handler {
		protected VCFHeader header;
		 boolean initialize(VCFHeader header) {
			this.header = header;
			return true;
			}
		abstract void visit(SortingCollection<Variation> sorter,VariantContext h) throws IOException;
		}

	private abstract class SOHandler extends Handler {
		@Override
		boolean initialize(VCFHeader h) {
			if(!super.initialize(h)) return false;
			
			return true;
			}
		@Override
		void visit(SortingCollection<Variation> sorter,VariantContext h)  throws IOException {
			
			}
		}

	
	@Override
	public int doWork(List<String> args) {
		try {
			final Path vcfpath = Paths.get(oneAndOnlyOneFile(args));
			IOUtil.assertFileIsReadable(vcfpath);
			
			
			try(VCFFileReader r= new VCFFileReader(vcfpath,true)) {
				final VCFHeader header= r.getFileHeader();
				List<GeneExtractor> extractors = new GeneExtractorFactory(header).getAllExtractors();
				try(CloseableIterator<VariantContext> iter = r.iterator()) {
					
					final SortingCollection<Variation> sorter = SortingCollection.newInstance(Variation.class, new RowCodec(),
							(A, B) -> A.compareVariantOnKeyAndPosition(B),
							writingSortingCollection.getMaxRecordsInRam(),
							writingSortingCollection.getTmpPaths());
					sorter.setDestructiveIteration(true);
					
					while(iter.hasNext()) {
						final VariantContext ctx= iter.next();
						for(GeneExtractor geneExtractor: extractors) {
							
							}
						
						
						}
					
					try (CloseableIterator<Variation> iter0 = sorter.iterator()) {
						try (EqualRangeIterator<Variation> iterE = new EqualRangeIterator<>(iter0, (A, B) -> A.compareVariantOnKey(B))) {
							while (iter.hasNext()) {
								final List<Variation> gene_variants = iterE.next();
								// sort and remove duplicates, if any
								Collections.sort(gene_variants,(A,B)->A.compareOnCoordinate(B));
								
								
							}
						}
					}
					
				
				}
			}
			
			
			return 0;
			}
		catch(Throwable err ) {
			LOG.error(err);
			return -1;
			}
		}

	public static void main(String[] args) {
		new RVTestsMakeSetFiles().instanceMainWithExit(args);
	}
}
