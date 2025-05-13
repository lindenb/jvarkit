/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.bgen;

import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

import com.github.lindenb.jvarkit.samtools.util.SimpleInterval;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.io.HtsPath;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.FeatureReader;
import htsjdk.tribble.TribbleIndexedFeatureReader;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;


public class BGenReader extends TribbleIndexedFeatureReader<BGenVariant,SeekableStream> {
	private static final Logger LOG = Logger.of(BGenReader.class).setDebug();
	private final SeekableStream seekableStream;
	//private final AbstractFeatureReader<BGenVariant, SeekableStream> featureReader;
	private final BGenCodec bgenCodec;
	/** auxiliary indexed VCF file that will be used as an coordinate index for the VCF */
	private Path filenameOrNull = null;
	private VCFFileReader vcfIndexFile = null;
	private final DelegateBGenReader delegate;
	
	private interface DelegateBGenReader extends Closeable {
		
		}
	
	private class  StreamBgenReader implements DelegateBGenReader {
		SeekableStream pbs;
		StreamBgenReader(InputStream in,final BGenCodec bgenCodec) {
			this.pbs = bgenCodec.makeSourceFromStream(in);
			}
		@Override
		public void close() throws IOException {
			
			}
		}
	
	private class AsbstractDelegateBGenReader extends AbstractFeatureReader<BGenVariant, SeekableStream> {
		AsbstractDelegateBGenReader(String path,final BGenCodec bgenCodec) {
			super(path,bgenCodec);
			}
		
		}

	
	public BGenReader(final InputStream in) throws IOException {
		this.bgenCodec = new BGenCodec(false);
		this.seekableStream = this.bgenCodec.makeSourceFromStream(in);
		this.bgenCodec.readHeader(this.seekableStream);
		}
	
	public BGenReader(final String path) throws IOException {
		this(SeekableStreamFactory.getInstance().getBufferedStream(SeekableStreamFactory.getInstance().getStreamFor(Objects.requireNonNull(path,"null input")),Defaults.NON_ZERO_BUFFER_SIZE));
		final HtsPath source = new HtsPath(path);
		if(source.isPath()) {
			this.filenameOrNull = source.toPath();
			}			
		}
	
	public BGenReader(final Path path) throws IOException {
		this(Objects.requireNonNull(path,"null input").toString());
		this.filenameOrNull = path;
		}
	
	public BGenHeader getHeader() {
		return this.bgenCodec.getHeader();
		}
	
	public Layout getLayout() {
		return getHeader().getLayout();
		}
	
	public Compression getCompression() {
		return getHeader().getCompression();
	}
	
	public long getSnpsOffset() {
		return this.bgenCodec.getSnpsOffset();
		}
	
	@Override
	public List<String> getSequenceNames() {
		return Collections.emptyList();
		}
	 
	/** change the reading pointer offset
	 * 
	 * @param position the physical offset
	 * @throws IOException the inutstream is not an instance of SeekableStream
	 */
	public void fseek(long position)  throws IOException {
		if(position<=0) throw new IllegalArgumentException("seek.pos<=0: "+position);
		this.seekableStream.seek(position);
		}
	

	
	private class Iter extends AbstractIterator<BGenVariant> implements
		CloseableTribbleIterator<BGenVariant>{
		private final boolean skipGenotypes;
		Iter(boolean skipGenotypes) {
			this.skipGenotypes=skipGenotypes;
			}
		
		@Override
		protected BGenVariant advance() {
			try {
				final BGenVariant v= BGenReader.this.readVariant();
				if(v==null) {
					return null;
					}
				else if(this.skipGenotypes) {
					BGenReader.this.skipGenotypes();
					return v;
					}
				else
					{
					return BGenReader.this.readGenotypes();
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public Iterator<BGenVariant> iterator() {
			return this;
			}
		@Override
		public void close() {
			
			}
		};
	
	public CloseableTribbleIterator<BGenVariant> iterator(boolean skipGenotypes) {
		return new Iter(skipGenotypes);
		}
	
	public CloseableTribbleIterator<BGenVariant> iterator() {
		return iterator(false);
		}
	
	private class QueryIterator extends AbstractIterator<BGenVariant> implements CloseableTribbleIterator<BGenVariant> {
		private Locatable query;
		private CloseableIterator<VariantContext> delegate;
		private boolean skipGenotypes;
		@Override
		protected BGenVariant advance() {
			try {
				for(;;) {
					if(delegate==null) return null;
					if(!delegate.hasNext()) {
						close();
						return null;
						}
					final VariantContext ctx = this.delegate.next();
					if(!ctx.overlaps(this.query)) continue;
					if(!ctx.hasAttribute(OFFSET_format_header_line.getID())) continue;
					final long offset= Long.parseLong(ctx.getAttributeAsString(OFFSET_format_header_line.getID(), "-1"));
					if(offset<=0L) throw new IllegalArgumentException("bad offset in "+ctx);
					BGenReader.this.fseek(offset);
					BGenVariant bv = BGenReader.this.readVariant();
					if(skipGenotypes) {
						BGenReader.this.skipGenotypes();
						}
					else
						{
						bv = BGenReader.this.readGenotypes();
						}
					return bv;
					}
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		@Override
		public void close()  {
			if(this.delegate!=null) this.delegate.close();
			this.delegate=null;
			}
		@Override
		public Iterator<BGenVariant> iterator() {
			return this;
			}
		}
	
	public CloseableTribbleIterator<BGenVariant> query(final String loc,int start,int end) {
		return query(new SimpleInterval(loc,start,end),false);
		}
	
	public CloseableTribbleIterator<BGenVariant> query(final Locatable loc, boolean skipGenotypes) {
		loadVcfIndex();
		final QueryIterator iter=new QueryIterator();
		iter.query=new SimpleInterval(loc);
		iter.delegate = Objects.requireNonNull(this.vcfIndexFile).query(loc);
		iter.skipGenotypes = skipGenotypes;
		return iter;
		}
	
	public boolean hasIndex() {
		if(this.vcfIndexFile!=null) return true;
		if(this.filenameOrNull==null)return false;
		final Path vcfPath = getVcfIndexName(this.filenameOrNull);
		if(!Files.exists(vcfPath)) return false;
		final String vcffilename = vcfPath.toAbsolutePath().toString();
		final  Path tbi = vcfPath.getFileSystem().getPath(ParsingUtils.appendToPath(vcffilename, FileExtensions.TABIX_INDEX));
		if(!Files.exists(tbi)) return false;
		try(VCFFileReader r=new VCFFileReader(vcfPath,false)) {
			final VCFHeader hder = r.getHeader();
			if(!hder.hasInfoLine(OFFSET_format_header_line.getID())) return false;
			}
		catch(Throwable err) {
			return false;
			}
		return true;
		}
	
	private void loadVcfIndex() {
		if(this.vcfIndexFile!=null) return;
		if(this.filenameOrNull==null) throw new RuntimeIOException("Cannot load a VCF index as the bgen is not a file (but probably a stream) ");
		final Path vcfPath = getVcfIndexName(this.filenameOrNull);
		if(!Files.exists(vcfPath))  throw new RuntimeIOException("cannot find associated VCF file "+vcfPath);
		this.vcfIndexFile=new VCFFileReader(vcfPath,true);
		final VCFHeader hder = this.vcfIndexFile.getHeader();
		if(!hder.hasInfoLine(OFFSET_format_header_line.getID())) {
			this.vcfIndexFile.close();
			this.vcfIndexFile=null;
			throw new RuntimeIOException("cannot use "+vcfPath+" as index because the header is missing ##INFO/"+OFFSET_format_header_line.getID());
			}
		}
	
	/** get the name of the VCF index for the bgen without checking it exists */
	public static Path getVcfIndexName(Path bgenPath) {
		final String filename = Objects.requireNonNull(bgenPath).toAbsolutePath().toString();
		final  String vcf = ParsingUtils.appendToPath(filename, VCF_INDEX_SUFFIX);
		return bgenPath.getFileSystem().getPath(vcf);
		}
	
	
	/** read the next variant or returns null if end of file */
	public BGenVariant readVariant() throws IOException {
		return this.bgenCodec.readVariant(this.seekableStream);
	}
	
	public void skipGenotypes() throws IOException{
		this.bgenCodec.skipGenotypes(this.seekableStream);
		}
	
	public BGenVariant readGenotypes()  throws IOException{
		return this.bgenCodec.readGenotypes(this.seekableStream);
	}
	
	
	
	/*
	private static List<int[]> buildPhasedIndexes(int n_alleles,int n_ploidy) {
		final List<int[]> container = new ArrayList<>(n_alleles*n_ploidy);
		buildPhasedIndexes(n_alleles,container,new int[n_ploidy],0);
		return container;
		}
	private static void buildPhasedIndexes(
			int n_alleles,
			List<int[]> container,
			int[] buffer,
			int ploidy_index
			) {
			if(ploidy_index==buffer.length) {
				container.add(Arrays.copyOf(buffer, buffer.length));
				}
			else
				{
				for(int a=0;a<n_alleles;++a) {
					int[] buffer2 =Arrays.copyOf(buffer, buffer.length);
					buffer2[ploidy_index]=a;
					buildPhasedIndexes(n_alleles,container,buffer2,ploidy_index+1);
					}
				}
			}*/
	
	@Override
	public void close()  {
		if(this.vcfIndexFile!=null) {
			vcfIndexFile.close();
			vcfIndexFile=null;
			}
		
		try {this.bgenCodec.close(this.seekableStream);}
		catch(Throwable err) {}
		}
	
	@Override
	public String toString() {
		return "BGenReader(file:"+this.filenameOrNull+")";
		}
}
