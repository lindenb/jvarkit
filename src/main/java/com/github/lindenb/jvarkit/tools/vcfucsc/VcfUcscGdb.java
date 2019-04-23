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

*/
package com.github.lindenb.jvarkit.tools.vcfucsc;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Path;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;


import org.apache.commons.jexl2.Expression;
import org.apache.commons.jexl2.JexlContext;
import org.apache.commons.jexl2.JexlException;
import org.apache.http.Header;
import org.apache.http.HttpEntity;
import org.apache.http.HttpStatus;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.client.methods.HttpHead;
import org.apache.http.impl.client.CloseableHttpClient;
import org.apache.http.impl.client.HttpClientBuilder;
import org.apache.http.impl.client.HttpClients;
import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BedFeature;
import org.broad.igv.bbfile.BigBedIterator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.broad.tribble.util.SeekableStream;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.LineIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import org.broad.tribble.util.SeekableFileStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFIterator;
/**
BEGIN_DOC

## input

chromosomes in the vcf must been in the same notation than the remote ucsc files (with a 'chr' prefix )

## Jexl expression

the following names are defined for the BED jexl context

 * **bed** : a **BedFeature** https://github.com/lindenb/bigwig/blob/master/src/org/broad/igv/bbfile/BedFeature.java
 * **length** : the length of the bigBed item
 * **contig** : the contig of the bigBed file
 * **chrom** : the contig of the bigBed file
 * **start** : the getStartBase of the bigBed file
 * **end** : the getEndBase of the bigBed file
 * **tokens** : a `List<String>` containing the columns of the bed line
 * **rest** : a `List<String>` containing the bed columns but the chrom/start/end
 * **line** : a `String` the whole bed line, separated with tab



END_DOC
*/
@Program(
		name="vcfucscgdb",
		description="annotate an VCF with ucsc remote bigbed/bigwig files.",
		keywords={"ucsc","mysql","bigwig","bigbed"},
		generate_doc=false,
		creationDate="20190423",
		modificationDate="20190423"
		)
public class VcfUcscGdb extends Launcher {
	private static final Logger LOG = Logger.build(VcfUcscGdb.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-r","--resource"},description=OPT_OUPUT_FILE_OR_STDOUT,required=true)
	private Path resourceFile = null;
	@Parameter(names={"-C"},description="try to convert chromosome name to match the UCSC nomenclature 1->'chr1'")
	private boolean forceUcscNames = false;
	@Parameter(names={"-f","--fraction"},description="fraction of the variant that should be covered by the bed feature [0-1] 1.0 : all variant bases must be covered.")
	private float fractionOverlap = 1.0f;

	
	
	private final HttpClientBuilder hb = HttpClients.custom();
	
	private static class NotIndexOutOfBoundList extends AbstractList<String> {
		final List<String> delegate;
		NotIndexOutOfBoundList(final List<String> delegate) {
			this.delegate= delegate;
			}
		@Override
		public String get(int index) {
			return index<0 || index >= this.delegate.size()?"":this.delegate.get(index);
			}
		@Override
		public int size() {
			return delegate.size();
			}
		@Override
		public String toString() {
			return String.join("\t", delegate);
			}
		}
	
	private abstract class LocatableJEXLContext implements JexlContext
		{
		@Override
		public void set(final String key, Object arg1) {
			throw new UnsupportedOperationException();
			}
		}
	
	private class WigItemJEXLContext extends LocatableJEXLContext
		{
		private final WigItem wigItem;
		WigItemJEXLContext(final WigItem wigItem) {
			this.wigItem = wigItem;
			}
		
		
		@Override
		public Object get(final String name) {
			if(name.equals("wig")) return this.wigItem;
			if(name.equals("length")) return CoordMath.getLength(this.wigItem.getStartBase(), this.wigItem.getEndBase());
			if(name.equals("contig")) return this.wigItem.getChromosome();
			if(name.equals("chrom")) return this.wigItem.getChromosome();
			if(name.equals("start")) return this.wigItem.getStartBase();
			if(name.equals("end")) return this.wigItem.getEndBase();
			if(name.equals("value")) return wigItem.getWigValue();
			LOG.warn("unknown property "+name);
			return null;
			}
		@Override
		public boolean has(final String key) {
			if(key.equals("wig")) return true;
			if(key.equals("length")) return true;
			if(key.equals("contig")) return true;
			if(key.equals("chrom")) return true;
			if(key.equals("start")) return true;
			if(key.equals("end")) return true;
			if(key.equals("value")) return true;
			return false;
			}
		
		@Override
		public String toString() {
			return String.join("|",
				this.wigItem.getChromosome(),
				String.valueOf(this.wigItem.getStartBase()),
				String.valueOf(this.wigItem.getEndBase()),
				String.valueOf(this.wigItem.getWigValue())
				);
			}
		}
	
	private class BedJEXLContext extends LocatableJEXLContext
		{
		private final BedFeature bedFature;
		private List<String> rest = null;
		private List<String> tokens = null;
		BedJEXLContext(final BedFeature bedFature) {
			this.bedFature = bedFature;
			}
		
		private List<String> getRest() {
			if(this.rest==null) {
				this.rest = new NotIndexOutOfBoundList(Arrays.asList(this.bedFature.getRestOfFields()));
				}
			return this.rest;
			}

		
		private List<String> getTokens() {
			if(this.tokens==null) {
				final List<String> R = getRest();
				final List<String> fields = new ArrayList<>(rest.size()+3);
				fields.add(this.bedFature.getChromosome());
				fields.add(String.valueOf(this.bedFature.getStartBase()));
				fields.add(String.valueOf(this.bedFature.getEndBase()));
				fields.addAll(R);
				this.tokens = new NotIndexOutOfBoundList(fields);
				}
			return this.tokens;
			}
		
		@Override
		public Object get(final String name) {
			if(name.equals("bed")) return this.bedFature;
			if(name.equals("length")) return CoordMath.getLength(this.bedFature.getStartBase(), this.bedFature.getEndBase());
			if(name.equals("contig")) return this.bedFature.getChromosome();
			if(name.equals("chrom")) return this.bedFature.getChromosome();
			if(name.equals("start")) return this.bedFature.getStartBase();
			if(name.equals("end")) return this.bedFature.getEndBase();
			if(name.equals("tokens")) return this.getTokens();
			if(name.equals("rest")) return getRest();
			if(name.equals("line")) return String.join("\t",this.getTokens());
			LOG.warn("unknown property "+name);
			return null;
			}
		@Override
		public boolean has(final String key) {
			if(key.equals("bed")) return true;
			if(key.equals("length")) return true;
			if(key.equals("contig")) return true;
			if(key.equals("chrom")) return true;
			if(key.equals("start")) return true;
			if(key.equals("end")) return true;
			if(key.equals("tokens")) return true;
			if(key.equals("rest")) return true;
			if(key.equals("line")) return true;
			return false;
			}
		@Override
		public String toString() {
			return String.join("|", getTokens());
			}
		}
	
	private abstract class AbstractJexlExprHandler {
		protected final Expression jexlExpr;
		private String exprStr;
		AbstractJexlExprHandler(final String exprStr) {
			this.exprStr = exprStr;
			this.jexlExpr = VariantContextUtils.engine.get().createExpression(exprStr);
			}
		protected Object eval(final LocatableJEXLContext t) {
			try {
				return this.jexlExpr.evaluate(t);
				}
			catch(final JexlException err) {
				throw new RuntimeException("Cannot evaluate JEXL expression \""+this.exprStr+"\" with BedRecord :"+t);
				}	
			}
	}
		
	
	private class JexlPredicate extends AbstractJexlExprHandler implements Predicate<LocatableJEXLContext> {
		JexlPredicate(String expr) {
			super(expr);
			}
		@Override
		public boolean test(final LocatableJEXLContext t) {
			final Object o = eval(t);
			if(o==null) throw new RuntimeException("JEXL expression \""+this.jexlExpr+"\" returned false");
			if(!(o instanceof Boolean)) throw new RuntimeException("JEXL expression \""+this.jexlExpr+"\" returned something that is not a boolean.");
			return Boolean.class.cast(o);
			}
		}
	
	private class JexlToString extends AbstractJexlExprHandler implements Function<LocatableJEXLContext,String> {
		JexlToString(String expr) {
			super(expr);
			}
		@Override
		public String apply(final LocatableJEXLContext t) {
			final Object o = eval(t);
			if(o==null) return "";
			return o.toString();
			}
		}
	
	
	
	private class RemoteBigFile
		implements Closeable
		{
		String name;
		String url;
		String description;
		int limit  = -1;
		private BBFileReader bbr =null;
		private SeekableStream seekableStream;
		private Predicate<LocatableJEXLContext> accept = L->true;
		private Function<LocatableJEXLContext,String> converter = L->L.toString();
		
		Set<String> query(String contig,int start,int end) throws IOException {
			final double variant_len = end-start+1;
			open();
			if(bbr.isBigBedFile()) {
				final Set<String>  annots = new LinkedHashSet<>();
				final BigBedIterator iter = bbr.getBigBedIterator(contig, start, contig, end, false /* contained */);
				while(iter.hasNext() && (this.limit==-1 || annots.size()< this.limit))
					{
					final BedFeature bed = iter.next();
					if(!bed.getChromosome().equals(contig)) continue;
					if(bed.getEndBase() < start) continue;
					if(bed.getStartBase() > end) continue;
					
					double overlap = CoordMath.getOverlap(bed.getStartBase(), bed.getEndBase(), start, end);
					
					if(overlap / variant_len <fractionOverlap) continue;
					
					final BedJEXLContext jexlCtx = new BedJEXLContext(bed);
					
					if(!accept.test(jexlCtx)) continue;
					final String annotation = this.converter.apply(jexlCtx);
					if(!StringUtils.isBlank(annotation)) annots.add(annotation);
					
					}
				return annots;
				}
			else if(bbr.isBigWigFile()) {
				final Set<String>  annots = new LinkedHashSet<>();
				final BigWigIterator iter = bbr.getBigWigIterator(contig, start, contig, end, false /* contained */);
				while(iter.hasNext() && (this.limit==-1 || annots.size()< this.limit))
					{
					final WigItem wiggle = iter.next();
					if(!wiggle.getChromosome().equals(contig)) continue;
					if(wiggle.getEndBase() < start) continue;
					if(wiggle.getStartBase() > end) continue;
					
					double overlap = CoordMath.getOverlap(wiggle.getStartBase(), wiggle.getEndBase(), start, end);
					
					if(overlap / variant_len <fractionOverlap) continue;
					
					final WigItemJEXLContext jexlCtx = new WigItemJEXLContext(wiggle);
					
					if(!accept.test(jexlCtx)) continue;
					final String annotation = this.converter.apply(jexlCtx);
					if(!StringUtils.isBlank(annotation)) annots.add(annotation);
					
					}
				return annots;
				}
			else
				{
				throw new IllegalStateException("not handled");
				}
			}
		
		void open() throws IOException {
			if(this.bbr == null) {
				LOG.info("opening "+this.url);
				if(IOUtil.isUrl(this.url))
					{
					final URL u = new URL(this.url);
					this.seekableStream = openRemoteSeekableStream(u);
					}
				else
					{
					final File file = new File(this.url);
					IOUtil.assertFileIsReadable(file);
					this.seekableStream = new SeekableFileStream(file);
					}
				
				this.bbr = new BBFileReader(url, seekableStream);
				}
			}
		
		@Override
		public void close() throws IOException {
			CloserUtil.close(this.seekableStream);
			CloserUtil.close(this.bbr);
			this.seekableStream=null;
			this.bbr=null;
			}
		}
	
private SeekableStream openRemoteSeekableStream(final URL url) 
    	throws IOException {
		
		final CloseableHttpClient httpClient = hb.build();
		
		final HttpHead httpHead = new HttpHead(url.toExternalForm());
		try
			{
			final CloseableHttpResponse response = httpClient.execute(httpHead);
			if (response.getStatusLine().getStatusCode() != HttpStatus.SC_OK ) {
				final String msg = "Unexpected Http status code "
     		            + response.getStatusLine()+" for "+ url
     		            ;
				response.close();
				httpClient.close();
 		        throw new IOException(msg);
 		   		}
			
			final Header contentLengthHeader = response.getFirstHeader("Content-Length");
			if(contentLengthHeader==null)
				{
				final String msg = "Cannot get Content-length for "+ url;
				response.close();
				httpClient.close();
 		        throw new IOException(msg);
				}
			
			long contentLength;
			try {
				contentLength = Long.parseLong(contentLengthHeader.getValue());
				if(contentLength<0) throw new NumberFormatException("Negative content length for "+url);
				}
			catch(final NumberFormatException err)
				{
				final String msg = "bad Content-Length in "+contentLengthHeader+" for "+ url;
				response.close();
				httpClient.close();
 		        throw new IOException(msg,err);
				}
			response.close();
			return new CustomRemoteSeekeableStream(httpClient,url,contentLength);
			}
		catch(final IOException err)
			{
			CloserUtil.close(httpClient);
			throw err;
			}	
    	}


    /** custom SeekableStream reading remote files */    
	private static class CustomRemoteSeekeableStream extends SeekableStream
		{
		private final CloseableHttpClient httpClient;
	    private long position = 0L;
	    private final URL url;
	    private final long contentLength;

	    CustomRemoteSeekeableStream(final CloseableHttpClient httpClient,final URL url,final long contentLength) {
	    	this.httpClient = httpClient;
	    	this.url = url;
	    	this.contentLength = contentLength;
	    	}
	    
		@Override
		public boolean eof() throws IOException {
	        return this.contentLength > 0 && this.position >= this.contentLength;
		}

		@Override
		public long length() {
			return contentLength;
		}

		@Override
		public long position() throws IOException {
	        return position;
		}

		@Override
		public void seek(long position) throws IOException {
	        this.position = position;
		}

		@Override
		public int read(byte[] buffer, int offset, int len) throws IOException {
	        if (offset < 0 || len < 0 || (offset + len) > buffer.length) {
	            throw new IndexOutOfBoundsException("Offset="+offset+",len="+len+",buflen="+buffer.length);
	        }
	        if (len == 0 ) {
	            return 0;
	        }
	        if (this.position == this.contentLength) {
	            return -1;
	        }

	        CloseableHttpResponse httpResponse = null;
	        InputStream is = null;
	        int n = 0;
	        try {
	        	final HttpGet httpGet = new HttpGet(this.url.toExternalForm());
	        	
	        	

	            long endRange = position + len - 1;
	            // IF we know the total content length, limit the end range to that.
	            if (contentLength > 0) {
	                endRange = Math.min(endRange, contentLength);
	            }
	            final String byteRange = "bytes=" + position + "-" + endRange;
	            
	            httpGet.addHeader("Range", byteRange);
	          
	            httpResponse = this.httpClient.execute(httpGet);
	            if (httpResponse.getStatusLine().getStatusCode() != HttpStatus.SC_PARTIAL_CONTENT ) {
					final String msg = "Unexpected Http status code "
	     		            + httpResponse.getStatusLine()+" for "+ url +" in range "+byteRange;
					httpResponse.close();
					httpResponse = null;
	 		        throw new IOException(msg);
	 		   		}
	            final HttpEntity entity = httpResponse.getEntity();
	            
	            is = entity.getContent();

	            while (n < len) {
	                int count = is.read(buffer, offset + n, len - n);
	                if (count < 0) {
	                    if (n == 0) {
	                        return -1;
	                    } else {
	                        break;
	                    }
	                }
	                n += count;
	            	}
	            this.position += n;

	            return n;
	        	}
	        finally {
	            CloserUtil.close(is);
	            CloserUtil.close(httpResponse);
	        }
	    }
		
		@Override
		public int read() throws IOException {
			final byte []tmp=new byte[1];
	    	read(tmp,0,1);
	    	return (int) tmp[0] & 0xFF; 
	    	}
		}

	private List<RemoteBigFile> readRemoteResources(final Path path) throws IOException
		{
		final List<RemoteBigFile> remoteBigFiles = new ArrayList<>();
		IOUtil.assertFileIsReadable(path);
		try(BufferedReader br  = IOUtil.openFileForBufferedReading(path)) {
			final HashMap<String, String> hash = new HashMap<>();
			final Function<String, String> required = (K)->{
				if(!hash.containsKey(K)) throw new RuntimeIOException("Key \""+K+"\" missing. Found: "+hash.keySet());
				final String v= hash.get(K).trim();
				if(StringUtils.isBlank(v)) throw new RuntimeIOException("Key \""+K+"\" is empty");
				return v;
				};
			try( LineIterator iter = new LineIterator(br)) {
				for(;;) {
					final String line = (iter.hasNext()?iter.next():null);
					
					
					if(StringUtils.isBlank(line)) {
						if(!hash.isEmpty()) {
							final RemoteBigFile bf= new RemoteBigFile();
							bf.name = required.apply("name");
							bf.url = required.apply("url");
							if(remoteBigFiles.stream()
								.anyMatch(R->R.name.equals(bf.name)))
								{
								bf.close();
								throw new RuntimeIOException("Duplicate remote resource: "+hash);
								}
							if(hash.containsKey("accept")) {
								bf.accept = new JexlPredicate(hash.get("accept"));
								}
							if(hash.containsKey("tostring")) {
								bf.converter = new JexlToString(hash.get("tostring"));
								}
							if(hash.containsKey("desc")) {
								bf.description = hash.get("desc");
								}
							else
								{
								bf.description = "Data from "+bf.url; 
								}
							if(hash.containsKey("limit")) {
								bf.limit  = Integer.parseInt(hash.get("limit"));
								}
							remoteBigFiles.add(bf);
							}
						if(line==null) break;
						hash.clear();
						continue;
						}
					if(line.startsWith("#")) continue;
					int sep = line.indexOf(':');
					if(sep==-1) sep =  line.indexOf('=');
					if(sep==-1) throw new RuntimeIOException("Cannot find ':' or '=' in  "+line);
					final String key = line.substring(0,sep).toLowerCase().trim();
					if(hash.containsKey(key)) throw new RuntimeIOException("Duplicate key "+key+" in resource: "+hash);
					final String value = line.substring(sep+1).trim();
					hash.put(key, value);
					}
				}
			}
		return remoteBigFiles;
		}
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iter, VariantContextWriter out) {
		List<RemoteBigFile> remoteBigFiles = new ArrayList<>();
		try {
			final ContigNameConverter contigNameConverter = 
					this.forceUcscNames
					? ContigNameConverter.createConvertToUcsc()
					: ContigNameConverter.getIdentity();

					remoteBigFiles = this.readRemoteResources(this.resourceFile);
			final VCFHeader header0 = iter.getHeader();
			
			final VCFHeader header = new VCFHeader(header0);
			for(final RemoteBigFile rsrc: remoteBigFiles) {
				final VCFInfoHeaderLine info = new VCFInfoHeaderLine(rsrc.name, VCFHeaderLineCount.UNBOUNDED,VCFHeaderLineType.String, rsrc.description);
				if(header.getInfoHeaderLine(rsrc.name)!=null) {
					LOG.error("Duplicate / already existing INFO/"+rsrc.name);
					return -1;
					}
				header.addMetaDataLine(info);
				}
			
			
			out.writeHeader(header);
			while(iter.hasNext()) {
				final VariantContext ctx = iter.next();
				
				
				
				final String ctg = contigNameConverter.apply(ctx.getContig());
				if(StringUtils.isBlank(ctg)) {
					out.add(ctx);
					continue;
					}
				
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
				for(final RemoteBigFile rsrc: remoteBigFiles) {
					vcb.rmAttribute(rsrc.name);
					final Set<String> annots = rsrc.query(ctg,ctx.getStart(),ctx.getEnd());
					if(!annots.isEmpty()) {
						vcb.attribute(rsrc.name, new ArrayList<>(annots));
						}
					}				
				out.add(vcb.make());
				}
			
			return -1;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			remoteBigFiles.forEach(R->CloserUtil.close(R));
			}
		}
	
	@Override
	public int doWork(final List<String> args) {
		if(fractionOverlap<=0 || fractionOverlap>1) {
			LOG.error("bad fraction overlap");
			return -1;
			}
		return doVcfToVcf(args, this.outputFile);
		}
	public static void main(final String[] args) {
		new VcfUcscGdb().instanceMain(args);
	}
}
