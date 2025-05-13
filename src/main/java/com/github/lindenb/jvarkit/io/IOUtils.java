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
package com.github.lindenb.jvarkit.io;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.PushbackInputStream;
import java.io.Reader;
import java.io.StringWriter;
import java.io.Writer;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.InflaterInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.iterator.LineIterators;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.LineReader;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.samtools.Defaults;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

public class IOUtils {
	/*
	private static abstract class AbstractErrorChecker implements BooleanSupplier
		{
		private final long intervalMillis=10*1000;//10 sec
		private long last = System.currentTimeMillis();
		protected abstract boolean checkNoError();
		private boolean is_ok=true;
		
		@Override
		public boolean getAsBoolean() {
			final long now = System.currentTimeMillis();
			if(is_ok && last+intervalMillis>now) {
				last=now;
				is_ok=checkNoError();
			}
			
			return is_ok;
			}
		
		}
	private static class PrintStreamChecker extends AbstractErrorChecker {
		private final PrintStream p;
		PrintStreamChecker(final PrintStream p) {
			this.p=p;
		}
		@Override
		protected boolean checkNoError() {
			return !p.checkError();
			}
		}
	private static class PrintWriterChecker extends AbstractErrorChecker {
		private final PrintWriter p;
		PrintWriterChecker(final PrintWriter p) {
			this.p=p;
		}
		@Override
		protected boolean checkNoError() {
			return !p.checkError();
			}
		}*/
	
	/** check directory is readeable and return the directory */
	public static Path assertDirectoryIsReadable(final Path dir) {
		IOUtil.assertDirectoryIsReadable(dir==null?null:dir.toFile());
		return dir;
	}
	
    /** Returns a default tmp directory. */
    public static File getDefaultTmpDir() {
    	return new File(System.getProperty("java.io.tmpdir","."));
    	}
    /** Returns a default tmp directory as a nio Path. */
    public static Path getDefaultTempDir() {
    	return Paths.get(System.getProperty("java.io.tmpdir","."));
    	}
    
	public static void copyTo(final File f,final OutputStream fous) throws IOException
		{
		copyTo(f.toPath(),fous);
		}
	
	public static void copyTo(final Path path,final OutputStream fous) throws IOException
		{
		InputStream fin=null;
		try {
			fin =Files.newInputStream(path);
			copyTo(fin,fous);
			fous.flush();
		} finally {
			CloserUtil.close(fin);
			}
		}
	
	
	/** copy one file to another*/
	public static void copyTo(final File fin,final File fout) throws IOException
		{
		if(fin.equals(fout)) {
			throw new IllegalArgumentException("copyTo src=dest:"+fin);
			}
		OutputStream fous=null;
		try {
			fous = Files.newOutputStream(fout.toPath());
			copyTo(fin,fous);
		} finally {
			CloserUtil.close(fous);
			}
		}

	public static void copyTo(final Path path,final Writer fous) throws IOException
		{
		Reader fin=null;
		try {
			fin = Files.newBufferedReader(path);
			copyTo(fin,fous);
			fous.flush();
		} finally {
			CloserUtil.close(fin);
			}
		}

	
	public static void copyTo(final File f,final Writer fous) throws IOException
		{
		copyTo(f.toPath(),fous);
		}
	
	public static void copyTo(final InputStream in,final File f) throws IOException
		{
		copyTo(in,f.toPath());
		}
	
	/** copy input stream to path, no gzip detection is performed */
	public static void copyTo(final InputStream in,final Path path) throws IOException
		{
		final OutputStream fous= Files.newOutputStream(path);
		copyTo(in,fous);
		fous.flush();
		fous.close();
		}
	/** copy reader  to path, no gzip detection is performed */
	public static void copyTo(final Reader in,final Path path) throws IOException
		{
		final Writer fous= Files.newBufferedWriter(path);
		copyTo(in,fous);
		fous.flush();
		fous.close();
		}
	
	/** copy reader  to file, no gzip detection is performed */
	public static void copyTo(final Reader in,final File out) throws IOException
		{
		copyTo(in,out.toPath());
		}

	
	/** copy content of stream `in` to `out` */
	public static void copyTo(final InputStream in,final OutputStream out) throws IOException
		{
		final byte buffer[]=new byte[2048];
		int nRead;
		while((nRead=in.read(buffer))!=-1)
			{
			out.write(buffer,0,nRead);
			}
		out.flush();
		}

	/** copy nbytes content of stream `in` to `out'. Error if cannot write all nbytes */
	public static void copyTo(InputStream in,OutputStream output ,long nbytes) throws IOException {
		long remains = nbytes;
		if(nbytes<0) throw new IllegalArgumentException("negative nbytes "+nbytes);
        final byte[] buffer = new byte[Defaults.NON_ZERO_BUFFER_SIZE];
        while(remains>0 ) {
        	final int toRead= (int)Math.min((long)buffer.length, remains);
        	final int nRead = in.read(buffer, 0, toRead);
        	if(nRead==-1) throw new IOException("Got "+nRead+" bytes while I expected "+toRead+" bytes");
            output.write(buffer, 0, nRead);
            remains-=nRead;
        	}
		}

	
	public static void copyTo(final Reader in,final Writer out) throws IOException
		{
		final char buffer[]=new char[2048];
		int nRead;
		while((nRead=in.read(buffer))!=-1)
			{
			out.write(buffer,0,nRead);
			}
		out.flush();
		}
	
	/** copy content of 'in' into String */
	public static String copyToString(final Reader in) throws IOException
		{
		final StringWriter sw = new StringWriter();
		copyTo(in, sw);
		return sw.toString();
		}
	
	
	/** compressed a String with gzip */
	public static byte[] gzipString(final String s) {
		byte array[] = s.getBytes();
		try
			 {
			 final ByteArrayOutputStream obj=new ByteArrayOutputStream(array.length);
		     final GZIPOutputStream gzip = new GZIPOutputStream(obj);
		     gzip.write(array);
		     gzip.finish();
		     gzip.flush();
		     gzip.close();	
		     obj.close();
		     return obj.toByteArray();
			 }
		catch(final IOException err) {
			throw new RuntimeIOException(err);
		}
	}
	
	public static boolean isRemoteURI(final String uri)
		{	
		if(uri==null) return false;
		if(!IOUtil.isUrl(uri)) return false;
	    return  uri.startsWith("http://") ||
				uri.startsWith("https://") ||
				uri.startsWith("ftp://")
				;
		}
	private static InputStream tryBGZIP(final InputStream in) throws IOException
		{
		final byte buffer[]=new byte[ 
				BlockCompressedStreamConstants.GZIP_BLOCK_PREAMBLE.length  ];
		
		final PushbackInputStream push_back=new PushbackInputStream(in,buffer.length+10);
		int nReads=push_back.read(buffer);
		push_back.unread(buffer, 0, nReads);
		
		try
			{
			if( nReads>= buffer.length && 
				buffer[0]==BlockCompressedStreamConstants.GZIP_ID1 &&
				buffer[1]==(byte)BlockCompressedStreamConstants.GZIP_ID2 &&
				buffer[2]==BlockCompressedStreamConstants.GZIP_CM_DEFLATE &&
				buffer[3]==BlockCompressedStreamConstants.GZIP_FLG &&
				buffer[8]==BlockCompressedStreamConstants.GZIP_XFL
				)
				{
				return new BlockCompressedInputStream(push_back);
				}
			}
		catch(final Exception err)
			{
			//not bzip
			}
		return new GZIPInputStream(push_back);
		}
	
	/** open URL for reading, decompressed if url ends with gz of bz2 . if uri is null or '-' , open stdin */
	public static InputStream openURIForReading(String uri) throws IOException
		{
		if(uri==null || uri.equals("-")) return System.in;
		if(isRemoteURI(uri))
			{
			final URL url=new URL(uri);
			InputStream in=url.openStream();
			//do we have .... azdazpdoazkd.vcf.gz?param=1&param=2
			int question=uri.indexOf('?');
			if(question!=-1) uri=uri.substring(0, question);
			if(isCompressedExtention(uri)) {
				if(uri.endsWith(".gz"))
					{
					return tryBGZIP(in);
					}
				else if(uri.endsWith(".bz2"))
					{
					return new BZip2CompressorInputStream(in);
					}
				}
			return in;
			}
		if(uri.startsWith("file://"))
			{
			uri=uri.substring(7);
			}
		return openPathForReading(Paths.get(uri));
		}
	
	/** open uri for buffered reader , expect UTF-8 */
	public static BufferedReader openURIForBufferedReading(final String uri) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openURIForReading(uri), Charset.forName("UTF-8")));
		}

	public static Reader openFileForReader(final File file) throws IOException
		{
		return openFileForReader(file.toPath());
		}
	
	/** open file with reader, detect with extension if file is compressed */
	public static Reader openFileForReader(final Path file) throws IOException
		{
		IOUtil.assertFileIsReadable(file);
		if(isCompressedExtention(file.getFileName().toString()))
			{
			return new InputStreamReader(openPathForReading(file));
			}
		return Files.newBufferedReader(file);
		}

	
	public static InputStream openFileForReading(final File file) throws IOException
		{
		return openPathForReading(file.toPath());
		}
	
	public static InputStream openPathForReading(final Path path) throws IOException
		{
		IOUtil.assertFileIsReadable(path);
		InputStream in= Files.newInputStream(path);
		if(path.getFileName().toString().endsWith(".gz"))
			{
			in = tryBGZIP(in);
			}
		else if(path.getFileName().toString().endsWith(".bz2")) {
			in  = new BZip2CompressorInputStream(in);
			}
		return Objects.requireNonNull(in,"cannot open "+path);
		}
	
	public static BufferedReader openPathForBufferedReading(final Path path) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openPathForReading(path), Charset.forName("UTF-8")));
		}
	
	public static BufferedReader openFileForBufferedReading(final File file) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(openFileForReading(file), Charset.forName("UTF-8")));
		}

    public static BufferedWriter openFileForBufferedWriting(final File file)  throws IOException
    	{
        return new BufferedWriter(new OutputStreamWriter(openFileForWriting(file)), Defaults.BUFFER_SIZE);
    	}
   
    
	/** return true if the file has a compressed suffix '.bfz' or '.gz' or '.bz2' */
	public static final boolean isCompressedExtention(final String s) {
		return s!=null && StringUtils.endsWith(s, ".gz",".bgz",".bz2");
		}

    
    /** return true if path has an interpretable compressed suffix */
    public static boolean isCompressed(final Path out) {
    		final String suff= Objects.requireNonNull(out).getFileName().toString();
    		return	isCompressedExtention(suff);
    		}
    
    /** output path for writing. The following extensions
     * are interpretted : vcf.gz, .bgz , .gz, .bz2 
     * @param file
     * @return
     * @throws IOException
     */
    public static OutputStream openPathForWriting(final Path file) throws IOException
		{
    	if(file==null) throw new IllegalArgumentException("path is null");
    	if(isCompressed(file)) {
	    	final String base = file.getFileName().toString();
		    if (base.endsWith(".vcf.gz") || base.endsWith(".bgz"))
		    	{
		        return new BlockCompressedOutputStream(
		        		file,
		        		BlockCompressedOutputStream.getDefaultCompressionLevel(),
		        		BlockCompressedOutputStream.getDefaultDeflaterFactory()
		        		);
		    	}
		    else if (base.endsWith(".bz2"))
		    	{
		        return new BZip2CompressorOutputStream(Files.newOutputStream(file));
		    	}
		    else if (base.endsWith(".gz"))
		    	{
		        return new GZIPOutputStream(Files.newOutputStream(file),true);
		    	}
		    else
		    	{
		    	throw new IllegalStateException("bad suffix ?? "+file);
		    	}
	    	}
	    else
	    	{
	        return Files.newOutputStream(file);
	    	}         
		}
    
    public static OutputStream openFileForWriting(final File file) throws IOException
    	{
    	return openPathForWriting(file.toPath());
    	}
    
    /** open a printwriter, compress if it ends with *.gz  */
    public static PrintWriter openFileForPrintWriter(final File file) throws IOException
		{
	    return openPathForPrintWriter(file.toPath());
		}
    
    /** open a printwriter, compress if it ends with *.gz  */
    public static PrintWriter openPathForPrintWriter(final Path file) throws IOException
		{
	    if (isCompressed(file))
	    	{
	        return new PrintWriter( openPathForWriting(file));
	    	}
	    else
	    	{
	        return new PrintWriter(Files.newBufferedWriter(file));
	    	}         
		}
    
    public static LineReader openFileForLineReader(File file) throws IOException
		{
    	return new SynchronousLineReader(openFileForReading(file));
		}
    
    /** @return a LineIterators that should be closed with CloserUtils */
    public static LineIterator openPathForLineIterator(final Path file) throws IOException
  		{
  		return  LineIterators.of(file);
  		}

    
    /** @return a LineIterators that should be closed with CloserUtils */
    public static LineIterator openFileForLineIterator(final File file) throws IOException
  		{
  		return openPathForLineIterator(file.toPath());
  		}
    
    public static LineReader openStdinForLineReader() throws IOException
		{
		return  new SynchronousLineReader(System.in);
		}
    
    public static BufferedReader openStdinForBufferedReader() throws IOException
		{
		return openStreamForBufferedReader(System.in);
		}
    /** convert input stream for buffered reader */
    public static BufferedReader openStreamForBufferedReader(final InputStream in) throws IOException
		{
		return  new BufferedReader(new InputStreamReader(in));
		}

    /** @return a LineIterators that should be closed with CloserUtils */
    public static LineIterator openStreamForLineIterator(final InputStream in) throws IOException
  		{
  		return  new LineIteratorImpl(openStreamForLineReader(in));
  		}
    
    /** @return a LineIterators that should be closed with CloserUtils */
    public static LineIterator openStdinForLineIterator() throws IOException
  		{
  		return  openStreamForLineIterator(System.in);
  		}

    public static LineReader openStreamForLineReader(final InputStream in) throws IOException
		{
		return  new SynchronousLineReader(in);
		}

    
    public static LineReader openURIForLineReader(String uri) throws IOException
		{
		return  new SynchronousLineReader(openURIForReading(uri));
		}
    /** @return a LineIterators that should be closed with CloserUtils */
    public static LineIterator openURIForLineIterator(String uri) throws IOException
  		{
  		return  new LineIteratorImpl(openURIForLineReader(uri));
  		}
    
    /** read String from DataInputStream
     * motivation: readUTF can't print lines larger than USHORTMAX
     *  */
    public static String readString(final DataInputStream in) throws IOException
    	{
    	int llength=in.readInt();
    	if(llength==-1) return null;
    	byte a[]=new byte[llength];
    	in.readFully(a, 0, llength);
    	return new String(a);
    	}
    
    /** write String to DataOutputStream
     * motivation: DataInputStream.readUTF can't print lines larger than USHORTMAX
     *  */
    public static void writeString(final DataOutputStream os,String s) throws IOException
		{
		if(s==null)
			{
			os.writeInt(-1);
			}
		else
			{
			byte array[]=s.getBytes();
			os.writeInt(array.length);
			os.write(array);
			}
		}
    public static final String UNROLL_FILE_MESSAGE="if filename ends with '.list' it is interpreted as a list of file (one file per line)";
	
    /** unroll one file if it ends with '.list' , else its a singleton of file. If file is null, returns empty list*/
    public static List<File> unrollFile(final File file)
		{
		return unrollPath(file==null?null:file.toPath()).
				stream().
				map(F->F.toFile()).
				collect(Collectors.toList())
				;
		}
    /** unroll one file if it ends with '.list' , else its a singleton of file. If file is null, returns empty list*/
    public static List<Path> unrollPath(final Path file)
		{
		if(file==null) return Collections.emptyList();
		IOUtil.assertFileIsReadable(file);
		if(!file.getFileName().toString().endsWith(".list")) return Collections.singletonList(file);
		try {
			return Files.readAllLines(file).stream().
				filter(L->!(StringUtils.isBlank(L) || L.startsWith("#"))).
				map(L->Paths.get(L)).
				collect(Collectors.toList());
			} 
		catch (final IOException e) {
			throw new RuntimeIOException(e);
			}
		}
    
    /** return 'inputs' as set, if a filename ends with '*.list'
     * is is considered as a file:one file per line
     * @param inputs files
     * @return set of files
     */
	public static LinkedHashSet<String> unrollFiles(java.util.Collection<String> inputs)
		{
		final LinkedHashSet<String> vcfFiles= new LinkedHashSet<>(inputs.size()+1);
		for(final String file : inputs)
			{
			if(file.isEmpty())
				{
				//ignore
				}
			else if(IOUtil.isUrl(file))
				{
				vcfFiles.add(file);
				}
			else if(file.endsWith(".list"))
				{
				File f = new File(file);
				IOUtil.assertFileIsReadable(f);
				 for (final String s : IOUtil.readLines(f))
				 	{
					if (s.endsWith("#")) continue;
					if (s.trim().isEmpty()) continue;
					vcfFiles.add(s);
					}
				}
			else
				{
				vcfFiles.add(file);
				}
			
			}
		return vcfFiles;
		}
	
    /** return 'inputs' as set, if a filename ends with '*.list'
     * is is considered as a file:one file per line
     * @param inputs files
     * @return set of files
     */
	public static LinkedHashSet<File> unrollFileCollection(java.util.Collection<File> inputs)
		{
		LinkedHashSet<File> vcfFiles= new LinkedHashSet<>(inputs.size()+1);
		for(final File file : inputs)
			{
			if(file.getName().endsWith(".list"))
				{
				IOUtil.assertFileIsReadable(file);
				 for (final String s : IOUtil.readLines(file))
				 	{
					if (s.endsWith("#")) continue;
					if (s.trim().isEmpty()) continue;
					vcfFiles.add(new File(s));
					}
				}
			else
				{
				vcfFiles.add(file);
				}
			
			}
		return vcfFiles;
		}
	
    /** 
     * new version of unrollFiles in 2018, for common usage...
     * only one list
     * duplicate file are removed
     * return List of Files
     */
	public static List<File> unrollFiles2018(final java.util.List<String> args)
		{
		if(args.isEmpty()) return Collections.emptyList();
		final LinkedHashSet<File> fileset = new LinkedHashSet<>();
		if(args.size()==1 && args.get(0).endsWith(".list"))
			{
			final File listFile = new File(args.get(0));
			IOUtil.assertFileIsReadable(listFile);
			IOUtil.readLines(listFile).forEach(s->{
				if (s.endsWith("#")) return;
				if (StringUtil.isBlank(s)) return;
				fileset.add(new File(s));
				});
			}
		else
			{
			fileset.addAll(args.stream().
					map(S->new File(S)).
					collect(Collectors.toList()/* to list to keep oreder */));
			}
		return new ArrayList<>(fileset);
		}
	
	 /** 
     * new version of unrollFiles
     * only one list
     * duplicate file are removed
     * return List of Files
     */
	public static List<Path> unrollPaths(final java.util.List<String> args)
		{
		if(args.isEmpty()) return Collections.emptyList();
		final LinkedHashSet<Path> fileset = new LinkedHashSet<>();
		if(args.size()==1 && args.get(0).equals("-")) /* stdin */ {
			return Collections.emptyList();
			}
		else if(args.size()==1 && args.get(0).endsWith(".list"))
			{
			final File listFile = new File(args.get(0));
			IOUtil.assertFileIsReadable(listFile);
			IOUtil.readLines(listFile).forEach(s->{
				if (s.endsWith("#")) return;
				if (StringUtil.isBlank(s)) return;
				fileset.add(Paths.get(s));
				});
			}
		else
			{
			fileset.addAll(args.stream().
					map(S->Paths.get(S)).
					collect(Collectors.toList()/* to list to keep oreder */));
			}
		return new ArrayList<>(fileset);
		}
	
	@Deprecated
	public static List<String> unrollStrings2018(final java.util.List<String> args) {
		return 	unrollStrings(args);
		}
	
	 /** 
     * if there is one argument and it ends with '.list' it is interpretted as a file containing a list of string
     * otherwise return a copy of 'args'
     */
	public static List<String> unrollStrings(final java.util.List<String> args)
		{
		if(args.size()==1 && args.get(0).equals("-")) {
			return new ArrayList<>();
			}
		else if(args.size()==1 && args.get(0).endsWith(".list"))
			{
			final Path listFile = Paths.get(args.get(0));
			IOUtil.assertFileIsReadable(listFile);
			try {
				return Files.lines(listFile).collect(Collectors.toList());
			} catch (IOException err) {
				throw new RuntimeIOException("Cannot unroll the content of file "+listFile,err);
				}
			}
		else
			{
			return new ArrayList<>(args);
			}
		}

	/** test wether the two first bytes are gzip */
	public static boolean isGZipCompressed(final byte[] twoBytes) {
		if ((twoBytes == null) || (twoBytes.length < 2)) {
			return false;
		} else {
			return    ((twoBytes[0] == (byte) (GZIPInputStream.GZIP_MAGIC))
					&& (twoBytes[1] == (byte) (GZIPInputStream.GZIP_MAGIC >> 8)));
		}
	}
	
	//http://stackoverflow.com/a/4818946/58082
	/** uncompress IF needed if the input stream is a gzipped stream 
	 * 
	 * */
	public static InputStream uncompress(final InputStream input) throws IOException {
		 if(input==null) throw new NullPointerException("input is null");
		 if(input instanceof GZIPInputStream) return input;
		 if(input instanceof BlockCompressedInputStream) return input;
		 final PushbackInputStream pb = new PushbackInputStream( input, 2 ); //we need a pushbackstream to look ahead
	     final byte [] signature = new byte[2];
	     final int nread=pb.read( signature ); //read the signature
	     if(nread!=-1) pb.unread( signature,0,nread); //push back the signature to the stream
	     if(nread==signature.length && isGZipCompressed(signature)) //check if matches standard gzip magic number
	       return new GZIPInputStream( pb );
	     else 
	       return pb;
	}
	
	
	/** converts a BufferedReader to a line Iterator. use return LineIterators.of(r);  */
	@Deprecated
	public static  LineIterator toLineIterator(final BufferedReader r)
		{
		return LineIterators.of(r);
		}
	
	/** Prevent an output stream to be closed. 
	 *  Wrap into a stream that will be flushed instead of being closed 
	 *  For example calling vcfwriter.close in a zip entry. */
	public static OutputStream uncloseableOutputStream(final OutputStream os) {
		return new FilterOutputStream(os) {
			@Override
			public void close() throws IOException {
				this.flush();
				os.flush();
				}
		};
	}
	
	/** safe flusher */
	public static void flush(final Writer w) {
		if(w==null) return;
		try { w.flush();} catch(Throwable err) {}
	}
	/** safe flusher */
	public static void flush(final OutputStream w) {
		if(w==null) return;
		try { w.flush();} catch(Throwable err) {}
	}
	
	/** write something in a file */
	public static void cat(final Object o,final Path out, boolean append)
		{
		PrintWriter w=null;
		try
			{
			w= new PrintWriter(
					append?
					Files.newBufferedWriter(out, StandardOpenOption.APPEND):
					Files.newBufferedWriter(out)
					);
			w.print(o);
			w.flush();
			w.close();
			w=null;
			}
		catch(final IOException err)
			{
			throw new RuntimeIOException(err);
			}
		finally
			{
			if(w!=null) w.close();
			}
		}

	
	/** write something in a file */
	public static void cat(final Object o,final File out, boolean append)
		{
		cat(o,out.toPath(),append);
		}
	
	/** create tmp directory into existing another directory */
	 public static File createTempDir(final String prefix, final String suffix,final File parentDir) {
	     IOUtil.assertDirectoryIsWritable(parentDir);   
		 try {
	            final File tmp = File.createTempFile(prefix, suffix,parentDir);
	            if (!tmp.delete()) {
	                throw new SAMException("Could not delete temporary file " + tmp);
	            }
	            if (!tmp.mkdir()) {
	                throw new SAMException("Could not create temporary directory " + tmp);
	            }
	            return tmp;
	        } catch (final IOException e) {
	            throw new RuntimeIOException("Exception creating temporary directory in "+parentDir, e);
	        }
	    }
	 
	 /** return the string after the last dot, INCLUDING THE DOT */
	public static final String getFileSuffix(final Path path) {
		if(path==null) throw new IllegalArgumentException("path is null");
		final String s = path.getFileName().toString();
		final int dot = s.lastIndexOf('.');
		if(dot==-1) throw new IllegalArgumentException("cannot find dot character in "+path);
		return s.substring(dot);
		}
	
    /** inflate==gonfler uncompress byte array that was compressed using the Deflater algorithm */
	public static byte[] inflate(byte compressedBytes[]) {
		try {
			final InflaterInputStream zipIn = new InflaterInputStream(new ByteArrayInputStream(compressedBytes));
			final ByteArrayOutputStream baos = new ByteArrayOutputStream();
			copyTo(zipIn, baos); 
			baos.close();  
			/* uncompressedBytes */ 
			return baos.toByteArray();
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	/** read all path content into string */
	public static String slurpPath(final Path p) {
		IOUtil.assertFileIsReadable(p);
		try {
			return new String(Files.readAllBytes(p));
		} catch (final IOException e) {
			throw new RuntimeIOException("Cannot read content of "+p, e);
			}
		}
	
	/** read all input stream into String*/
	public static String slurpInputStream(final InputStream in) {
		try {
			try(Reader r=new InputStreamReader(in)) {
				return slurpReader(r);
			}
		} catch (final IOException e) {
			throw new RuntimeIOException("Cannot read content of inputStream", e);
			}
		}
	
	/** read all input stream into String*/
	public static String slurpReader(final Reader r) {
		try {
			return copyToString(r);
		} catch (final IOException e) {
			throw new RuntimeIOException("Cannot read content of Reader", e);
			}
		}

	
	 /** Returns the name of the file minus the common suffixes */
    public static String getFilenameWithoutCommonSuffixes(final Path f) {
        String s = f.getFileName().toString();
        final String suffixes[]= {
        		".bed",".txt",".zip",".bam",".cram",".sam",".vcf",".bcf",
        		".xls",".tsv",".csv",".gff",".gtf",".csi",".tbi",
        		".fasta",".fa",
        		".fastq",".fq",
        		".bb",".bw",".bigbed",".bigwig",".bigWig",
        		".tar",".gz",
        		".pdf",".xml",
        		".interval_list",
        		".sql",//UCSC
        		".sortedByCoord.out"//STAR
        		};
       boolean changed=true;
       while(changed) {
    	   changed = false;
        	for(int i=0;i< suffixes.length;i++) {
        		final String suff= suffixes[i];
        		if(s.equals(suff)) continue;
        		if(!s.endsWith(suff)) continue;
        		s = s.substring(0,s.length()-suff.length());
        		changed=true;
        		}
        	}
		return s;
    	}
	
    /** Returns the name of the file minus the directory and the extension (i.e. text after the last "." in the filename). 
     * repeat the procedure until 'repeat'. repeat==-1 : forever. '.' at beggining of the line is always kept (hidden linux file ) */
    public static String getFilenameWithoutSuffix(final Path f,int repeat) {
        String s = f.getFileName().toString();
        if(repeat<-1) throw new IllegalArgumentException("repeat <-1");
        for(;;) {
			if(repeat==0) return s;
			final int index = s.lastIndexOf('.');
			if(index<1) return s;//0 is for hidden file
			s = s.substring(0,index);
			if(repeat!=-1) repeat--;
        	}
    	}
    /** return default user agent for http connnection */
    public static String getDefaultUserAgent() {
    	return "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:137.0) Gecko/20100101 Firefox/137.0";
    	}
    
    /** Reads len bytes in a loop.
     * @param in The InputStream to read from
     * @param buf The buffer to fill
     * @param off offset from the buffer
     * @param len the length of bytes to read
     * @throws IOException if it could not read requested number of bytes 
     * for any reason (including EOF)
     */
    public static void readFully(final InputStream in, final byte buf[], int off, int len ) throws IOException {
      int toRead = len;
      while ( toRead > 0 ) {
        int ret = in.read( buf, off, toRead );
        if ( ret < 0 ) {
          throw new IOException( "Premeture EOF from inputStream");
        }
        toRead -= ret;
        off += ret;
      }
    }
    /** Reads len bytes in a loop.
     * @param in The InputStream to read from
     * @param buf The buffer to fill
     * @throws IOException if it could not read requested number of bytes 
     * for any reason (including EOF)
     */
    public static void readFully(final InputStream in, final byte buf[]) throws IOException {
    	readFully(in, buf,0,buf.length);
    }
    
    /** Similar to readFully(). Skips bytes in a loop.
     * @param in The InputStream to skip bytes from
     * @param len number of bytes to skip.
     * @throws IOException if it could not skip requested number of bytes 
     * for any reason (including EOF)
     */
    public static void skipFully(final InputStream in, long len ) throws IOException {
      while ( len > 0 ) {
        long ret = in.skip( len );
        if ( ret < 0 ) {
          throw new IOException( "Premeture EOF from inputStream");
        }
        len -= ret;
      }
    }
    
    /**
     * call mayBeGzippedInputStream with recursive=true
     * @param in inputstream
     * @return
     */
    
    public static InputStream mayBeGzippedInputStream(InputStream in) throws IOException {
    	return mayBeGzippedInputStream(in,true);
    	}
    
    
    /**
     * uncompressed gzipped input stream if needed
     * @param in inputstream
     * @param recursive gunzip as long as the input stream is compressed
     * @return
     */
    
    public static InputStream mayBeGzippedInputStream(InputStream in,boolean recursive) throws IOException {
    	/* NO don't do this if gzipped twice ! */
	    	//if(in instanceof GZIPInputStream) return in;
	    	//if(in instanceof BlockCompressedInputStream) return in;
        for(;;) {
	    	// wrap the input stream into a BufferedInputStream to reset/read a BCFHeader or a GZIP
	        // buffer must be large enough to contain the BCF header and/or GZIP signature
	       final  BufferedInputStream  bufferedinput = new BufferedInputStream(in,1+ IOUtil.GZIP_HEADER_READ_LENGTH);
	        // test for gzipped inputstream
	        if(IOUtil.isGZIPInputStream(bufferedinput)) {
	            // this is a gzipped input stream, wrap it into GZIPInputStream
	            // and re-wrap it into BufferedInputStream so we can test for the BCF header
	            in = new GZIPInputStream(bufferedinput);
	            if(!recursive) return in;
	        	}
	        else
		        {
		        return bufferedinput;
		        }
	        }
    	}
	}
