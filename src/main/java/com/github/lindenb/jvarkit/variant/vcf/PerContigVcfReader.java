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
package com.github.lindenb.jvarkit.variant.vcf;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.variant.vcf.VCFReader;

/**
 * An instance of VCFReader initialized with a list of VCF where each vcf holds a specific chromosome
 */
public class PerContigVcfReader implements VCFReader {
	private final  VCFHeader firstHeader;
	private String currentContig = null;
	private ContigReader currentReader=null;
	
	/** information about one vcf file */
	private static class ContigReader  {
		final Path path;
		VCFReader vcfReader = null;
		final Set<String> contigs=new HashSet<>();
		ContigReader(final Path path) {
			this.path=path;
			}
		boolean hasContig(final String chrom) {
			return contigs.contains(chrom);
			}
		void open() {
			if(vcfReader!=null) return;
			vcfReader = new VCFFileReader(this.path, true);
			}
		public void close()  {
			if(vcfReader!=null) try{vcfReader.close();}catch(final Throwable err) {}
			vcfReader = null;
			}
		}
	/** all VCF files */
	private final List<ContigReader> contigReaders;
	/** map contig to VCF files */
	private final Map<String,ContigReader> contig2reader;

	/** 
	 * @param path a file with the .list extension containing the path to the indexed VCF
	 * @throws IOException
	 */
	public PerContigVcfReader(final Path path) throws IOException {
		if(!path.getFileName().toString().endsWith(".list")) {
			throw new IOException("filename "+path+" should end with .list");
			}
		final List<String> lines = Files.readAllLines(path);
		this.contigReaders = new ArrayList<>(lines.size());
		this.contig2reader = new HashMap<>(lines.size());
		for(String line: lines) {
			if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
			final Path vcfPath = Paths.get(line);
			IOUtil.assertFileIsReadable(vcfPath);
			final ContigReader ctgReader = new ContigReader(vcfPath);
			// use tabix to collect chromosomes
			try(TabixReader tabix = new TabixReader(line)) {
				ctgReader.contigs.addAll(tabix.getChromosomes());
				}
			// no variant is that vcf
			if(ctgReader.contigs.isEmpty()) continue;
			this.contigReaders.add(ctgReader);
			for(String ctg:ctgReader.contigs) {
				if(contig2reader.containsKey(ctg)) {
					throw new IOException("in "+path+" contig "+ctg+" is present in "+vcfPath+" and "+contig2reader.get(ctg));
					}
				contig2reader.put(ctg, ctgReader);
				}
			}
		if(this.contigReaders.isEmpty()) {
			throw new IOException("No vcf was found in "+path);
			}
		this.currentReader = this.contigReaders.get(0);
		this.currentReader.open();
		this.currentContig = this.currentReader.contigs.iterator().next();
		this.firstHeader = this.currentReader.vcfReader.getHeader();
		}
	
	
	@Override
	public void close() {
		for(ContigReader r:this.contigReaders) {
			r.close();
			}
		this.contigReaders.clear();
		this.contig2reader.clear();
		this.currentReader=null;
		this.currentContig=null;
		}

	@Override
	public VCFHeader getHeader() {
		return firstHeader;
	}

	@Override
	public CloseableIterator<VariantContext> query(final String chrom, int start, int end) {
		if(currentContig!=null && currentContig.equals(chrom)) {
			if(currentReader==null) //no vcf for this contig
				{
				return AbstractCloseableIterator.empty();
				}
			else
				{
				return currentReader.vcfReader.query(chrom, start, end);
				}
			}
		this.currentContig=chrom;
		// current opened VCF is not able to read chromosome
		if(currentReader!=null && !currentReader.hasContig(chrom)) {
			currentReader.close();
			}
		
		this.currentReader = contig2reader.get(chrom);
		if(currentReader!=null) {
			currentReader.open();
			}
		return this.query(chrom,start,end);//yes recursion
		}

	@Override
	public boolean isQueryable() {
		return true;
	}

	@Override
	public CloseableIterator<VariantContext> iterator() {
		final List<VCFIterator> iterators= new ArrayList<>(this.contigReaders.size());
		final VCFIteratorBuilder vib = new VCFIteratorBuilder();
		try {
			for(final ContigReader cr:this.contigReaders) {
				iterators.add(vib.open(cr.path));
				}
			}
		catch(IOException err) {
			for(VCFIterator t:iterators) {
				t.close();
				}
			throw new RuntimeIOException(err);
			}
		return new SequentialIterator(iterators);
		}

	private static class SequentialIterator extends AbstractCloseableIterator<VariantContext> {
		private final List<VCFIterator> iterators;
		private int idx=0;
		SequentialIterator(final List<VCFIterator> iterators) {
			this.iterators = iterators;
			}
		@Override
		protected VariantContext advance() {
			while(idx < this.iterators.size()) {
				if(this.iterators.get(idx).hasNext()) {
					return this.iterators.get(idx).next();
					}
				this.iterators.get(idx).close();
				++idx;
				}
			return null;
			}
		
		@Override
		public void close() {
			for(VCFIterator t:iterators) {
				t.close();
				}
			iterators.clear();
			idx=0;
			}
		}
	
}
