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
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.LongList;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.seekablestream.ISeekableStreamFactory;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.BinaryCodec;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

/** class containing method to access variants by index */
public class VCFByIndex implements LongList<VariantContext>, AutoCloseable {
	private static Logger LOG=Logger.of(VCFByIndex.class);
	public static final String INDEX_SUFFIX =".ith";
	private final AbstractSeekStream vcfStream;
	private final SeekableStream indexStream;
	private final long n_items;
	private final BinaryCodec binaryCodec;
	private final VCFUtils.CodecAndHeader codecHeader;
	
	private static  abstract class AbstractSeekStream extends InputStream{
		protected String path;
		private final StringBuilder buffer=new StringBuilder();
		AbstractSeekStream(final String path) {
			this.path = path;
			}
		public abstract long getPosition()  throws IOException ;
		public abstract void fseek(long position)throws IOException;
		@Override
		public final int read(byte[] b) throws IOException {
			return read(b,0,b.length);
			}
		
		/** read an return the next line, return null on EOF */
		String readLine() throws IOException
			{
			this.buffer.setLength(0);
			int c;
		    while ((c = this.read()) >= 0 && c != '\n')
		        buffer.append((char) c);
		    if (c < 0) return null;
		    return buffer.toString();
			}
		
		
		/** just read the line without building a string , return false when there is no more line */
		private  boolean nextLine() throws IOException {
			int c;
		    while ((c = this.read()) >= 0 && c != '\n') {
		    	}
		    return c>=0;
		    }
		
		
		/** read all the lines until the VCF CHROM line is found */
		VCFUtils.CodecAndHeader parseHeaderAndCodec() throws IOException {
			String line=null;
			final List<String> headerLines=new ArrayList<>();
			while((line=readLine())!=null)
				{
				headerLines.add(line);
				if(line.startsWith("#CHROM"))  break;
				}
			return VCFUtils.parseHeader(headerLines);
			}	
		
		@Override
		public String toString() {
			return path;
			}
		}
	
	private static class PlainVcfSeekStream extends AbstractSeekStream {
		private SeekableStream in;
		PlainVcfSeekStream(final ISeekableStreamFactory ssf,final String path) throws IOException {
			super(path);
			this.in = ssf.getBufferedStream(ssf.getStreamFor(path));
			}
		
		
		@Override
		public long getPosition() throws IOException {
			return this.in.position();
			}
		@Override
		public void fseek(long position) throws IOException {
			this.in.seek(position);
			}
		@Override
		public int read() throws IOException {
			return this.in.read();
			}
		@Override
		public int read(byte[] b, int off, int len) throws IOException {
			return this.in.read(b, off, len);
			}
		}
	private static class BGZVcfSeekStream extends AbstractSeekStream {
		final BlockCompressedInputStream in;
		BGZVcfSeekStream(final ISeekableStreamFactory ssf,final String path) throws IOException {
			super(path);
			final SeekableStream in0 = ssf.getBufferedStream(ssf.getStreamFor(path));
			this.in = new BlockCompressedInputStream(in0);
			}
		@Override
		public void fseek(long position) throws IOException {
			this.in.seek(position);
			}
		@Override
		public long getPosition() {
			return this.in.getPosition();
			}
		@Override
		public int read() throws IOException {
			return this.in.read();
			}
		@Override
		public int read(byte[] b, int off, int len) throws IOException {
			return in.read(b, off, len);
			}
		@Override
		public void close() throws IOException {
			in.close();
			}
		}
	
	public VCFByIndex(final Path vcfPath) throws IOException {
		this(vcfPath.toString());
		}
	public VCFByIndex(final Path vcfPath,final Path idxPath) throws IOException {
		this(vcfPath.toString(),idxPath.toString());
		}

	public VCFByIndex(final String vcfPath) throws IOException {
		this(vcfPath,ParsingUtils.appendToPath(vcfPath, INDEX_SUFFIX));
		}
	
	public VCFByIndex(final String vcfPath,final String indexPath) throws IOException {
		final ISeekableStreamFactory ssf= SeekableStreamFactory.getInstance();		
		this.vcfStream = openVcfFile(ssf,vcfPath);
		this.codecHeader = this.vcfStream.parseHeaderAndCodec();

		this.indexStream = ssf.getBufferedStream(ssf.getStreamFor(indexPath));
		this.n_items = this.indexStream.length()/Long.BYTES;
		this.binaryCodec = new BinaryCodec(this.indexStream);
		}
	
	public VCFHeader getHeader()  {
		return this.codecHeader.header;
		}
	
	@Override
	public VariantContext get(long idx) {
		if(idx<0 || idx>=this.size()) throw new IndexOutOfBoundsException("0<="+idx+"<"+size());
		try {
			this.indexStream.seek(idx*Long.BYTES);
			final long offset = this.binaryCodec.readLong();
			this.vcfStream.fseek(offset);
			final String line = this.vcfStream.readLine();
			if(line==null) throw new NullPointerException("cannot readline at offset="+offset);
			return this.codecHeader.codec.decode(line);
			}
		catch(IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	
	@Override
	public long size() {
		return this.n_items;
		}
	
	@Override
	public void close() {
		try { this.vcfStream.close();} catch(Throwable err) {}
		try { this.binaryCodec.close();} catch(Throwable err) {}
		try { this.indexStream.close();} catch(Throwable err) {}
		}

	@Override
	public String toString() {
		return "VCFByIndex("+this.vcfStream.path+")";
		}
	
	/** return the index path for the vcf without creating it */
	public static Path getIndexForPath(final Path vcf) {
		return Paths.get(vcf.toString()+INDEX_SUFFIX);
		}
	
	public static Path buildIndex(final Path path) throws IOException {
		return buildIndex(path, getIndexForPath(path));
		}
	public static Path buildIndex(Path vcfpath, final Path indexPath) throws IOException {
		IOUtil.assertFileIsReadable(vcfpath);
		if(Files.isSameFile(vcfpath, indexPath)) throw new IllegalArgumentException("same file "+vcfpath);
		return buildIndex(vcfpath.toAbsolutePath().toString(), indexPath);
		}
	
	/** return true if index can be opened */
	public static boolean indexExists(final String f) {
		try(SeekableStream st=SeekableStreamFactory.getInstance().getStreamFor(f)) {
			st.eof();//remove warning st no used
			return true;
			}
		catch(IOException err) {
			return false;
			}
		}
	
	/** build index from path to file indexPath, return the indexPath */
	public static Path buildIndex(final String vcfpath, final Path indexPath) throws IOException {
		if(vcfpath.endsWith(FileExtensions.BCF)) {
			throw new IllegalArgumentException("cannot build a BCF from " + vcfpath);
			}
		if(!indexPath.getFileName().toString().endsWith(INDEX_SUFFIX)) {
			LOG.warning(indexPath.toString()+": suffix is not "+INDEX_SUFFIX);
			}
		final Path tmpIndex = Files.createTempFile("tmp.", INDEX_SUFFIX);
	
		try(AbstractSeekStream in = openVcfFile(SeekableStreamFactory.getInstance(),vcfpath)) {
			/* parse header */
			in.parseHeaderAndCodec();
			
			try (BinaryCodec codec=new BinaryCodec(tmpIndex, true)) {
				long currentPos = in.getPosition();
				/** loop over each line */
				while(in.nextLine()) {
					codec.writeLong(currentPos);
					currentPos = in.getPosition();
					}
				codec.getOutputStream().flush();
				}
			
			IOUtils.copyTo(tmpIndex, indexPath);
			return indexPath;
		}
		finally {
			Files.deleteIfExists(tmpIndex);
		}
	}
	
	
	private static AbstractSeekStream openVcfFile(final ISeekableStreamFactory ssf,final String vcfFile) throws IOException {
		if(IOUtil.hasBlockCompressedExtension(vcfFile))
			{
			return new BGZVcfSeekStream(ssf,vcfFile);
			}
		else if(vcfFile.endsWith(".vcf"))
			{
			return new PlainVcfSeekStream(ssf,vcfFile);
			}
		else
			{
			throw new IOException("Not a .vcf or .vcf.gz, .vcf.bgz file: "+vcfFile);
			}
		}
	

		

	
	
	
	
		
	



}
