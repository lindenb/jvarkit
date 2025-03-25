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
package com.github.lindenb.jvarkit.tools.jbrowse2;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.io.Writer;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.net.ContentType;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;

import htsjdk.samtools.SamFiles;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC


END_DOC
*/
public class JBrowse2Server  extends Launcher {
	private static final Logger LOG = Logger.build(JBrowse2Server.class).make();
	@Parameter(names="--zip",description="JBrowse2 archive source")
	private String jbrowse2url = "https://github.com/GMOD/jbrowse-components/releases/download/v3.0.3/jbrowse-web-v3.0.3.zip";
	@Parameter(names="--port",description="server port.")
	private int serverPort = 8080;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required = true)
	private Path reference;

	private static abstract class AbstractContent {
		protected void write(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
			}
		}
	private static class JsonConfigFile extends AbstractContent {
		final JsonObject config;
		JsonConfigFile(final JsonObject config) {
			this.config = config;
			}
		@Override
		protected void write(HttpServletRequest request, HttpServletResponse response)
				throws ServletException, IOException {
			response.setContentType(ContentType.APPLICATION_JSON.getMimeType());
			final Writer w=response.getWriter();
			new Gson().toJson(this.config, w);
			w.flush();
			}
		}
	private static class VirtualFile extends AbstractContent {
		final Path path;
		final long file_length;
		VirtualFile(final Path path) throws IOException {
			this.path=path;
			this.file_length = Files.size(path);
			}
		
		@Override
		protected void write(HttpServletRequest request, HttpServletResponse response)
				throws ServletException, IOException {
			response.setContentType(ContentType.APPLICATION_OCTECT_STREAM.getMimeType());
			String range=request.getHeader("Range");
			if(!StringUtils.isBlank(range) ) {
				if(!range.startsWith("bytes="))throw new IOException("range ?"+range);
				final long[] offset=Arrays.stream(CharSplitter.HYPHEN.split(range.substring(6))).
						mapToLong(Long::valueOf).
						map(v->Math.min(v,file_length)).
						toArray();
				if(offset.length!=2) throw new IOException("range ?"+range);
				final long out_len = 1+offset[1]-offset[0];
				response.setContentLengthLong(this.file_length);
				final OutputStream out= response.getOutputStream();
				
				System.err.println(range+" len="+out_len+" ss="+Arrays.toString(offset)+" length::"+this.file_length);
				
				try(RandomAccessFile rnd=new RandomAccessFile(this.path.toFile(),"r")) {
					long remain=out_len;
					rnd.seek(offset[0]);
					final byte buffer[]=new byte[2048];
					while(remain>0L) {
						final int nRead=rnd.read(buffer,0,(int)Math.min((long)buffer.length,remain));
						System.err.println("len:"+out_len+" remain:"+remain+" "+nRead);
						if(nRead==-1) break;
						out.write(buffer,0,nRead);
						remain -= nRead;
						}
					out.flush();
					}
				return;
				}
			
			response.setContentLengthLong(file_length);
			try(InputStream is=Files.newInputStream(this.path)) {
				IOUtils.copyTo(is, response.getOutputStream());
				}
			}
		}
	private static class FileInZip extends AbstractContent {
		final String name;
		final byte[] content;
		FileInZip(String name,byte[] content) {
			this.name=name;
			this.content = content;
			}
		protected void write(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
			final OutputStream out = response.getOutputStream();
			final ContentType contentType = ContentType.fromSuffix(name).orElse(null);
			if(contentType!=null) {
				response.setContentType(contentType.getMimeType());
				}
			else
				{
				LOG.info("undefined extension for "+name);
				}
			response.setContentLengthLong(this.content.length);
			out.write(this.content, 0, this.content.length);
			out.flush();
			}
		}
	
	private static class JBrowse2Servlet extends HttpServlet {	 
		private static final long serialVersionUID = 1L;
		private final Map<String,AbstractContent> page2content;
		
		JBrowse2Servlet(
				final Map<String,AbstractContent> page2content
				) {
			this.page2content= Collections.unmodifiableMap(page2content);
			}
		private void debug(HttpServletRequest request) {
			Enumeration<String> att= request.getAttributeNames();
			while(att.hasMoreElements()) {
				System.err.println(" att="+att.nextElement());
				}
			att= request.getHeaderNames();
			while(att.hasMoreElements()) {
				String s=att.nextElement();
				System.err.println(" header="+s+" = "+request.getHeader(s));
				}
			
			System.err.println("METHOD:"+request.getMethod());
			System.err.println("PATH INFO:"+request.getPathInfo());
			System.err.println("QUERY STRING:"+request.getQueryString());
			System.err.println("QUERY CTX PATH:"+request.getContextPath());
			System.err.println();
			if(!(request.getMethod().equalsIgnoreCase("get") || request.getMethod().equalsIgnoreCase("post"))) {
				System.exit(-1);
				}
			}
		
		@Override
		protected void doGet(HttpServletRequest req, HttpServletResponse resp) throws ServletException, IOException {
			doPost(req, resp);
			}
		@Override
		protected void doPost(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
			
			String pathInfo = request.getPathInfo();		
			if(pathInfo.endsWith("/") || StringUtils.isBlank(pathInfo)) pathInfo  = "/index.html";
			if(!(pathInfo.endsWith(".js") || pathInfo.endsWith(".map"))) debug(request);
			if(this.page2content.containsKey(pathInfo)) {
				final AbstractContent page = this.page2content.get(pathInfo);
				page.write(request,response);
				
				return;
				}
			
			throw new IOException("request "+request.getContextPath()+" "+pathInfo);
			}
		}

	private abstract class FileTypeHandler {
		protected String getBuildName() {
			return IOUtils.getFilenameWithoutCommonSuffixes(JBrowse2Server.this.reference); 
			}
		protected JsonObject createUriLocation(String path) {
			final JsonObject idm7 = new JsonObject();
			idm7.addProperty("uri",path);
			idm7.addProperty("locationType","UriLocation");
			return idm7;
			}
		
		abstract boolean handle(final String filename,final Map<String,AbstractContent> page2content,JsonArray tracks) throws IOException;
		}

	private  class BamTypeHandler extends FileTypeHandler {
		@Override
		boolean handle(String filename,final Map<String,AbstractContent> page2content,JsonArray tracks)  throws IOException {
			if(!filename.endsWith(FileExtensions.BAM)) return false;
			final SamReaderFactory srf = SamReaderFactory.makeDefault();
			srf.referenceSequence(JBrowse2Server.this.reference);
			try {
				final Path bamPath =Paths.get(filename);
				IOUtil.assertFileIsReadable(bamPath);
				try(SamReader sr=srf.open(bamPath)) {
					if(!sr.hasIndex()) {
						throw new IOException("index missing for "+bamPath);
						}
					String name0 = sr.getFileHeader().getReadGroups().stream().
							map(RG->RG.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().
							orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
					String bam_name=name0+FileExtensions.BAM;
					String bai_name=name0+FileExtensions.BAM+FileExtensions.BAI_INDEX;
					int n=1;
					for(;;) {
						if(!page2content.containsKey(bam_name) && !page2content.containsKey( bai_name)) break;
						n++;
						bam_name=name0+"."+n+FileExtensions.BAM;
						bai_name=name0+"."+n+FileExtensions.BAM+FileExtensions.BAI_INDEX;
						}
					VirtualFile virtBam = new VirtualFile(bamPath);
					page2content.put("/"+bam_name, virtBam);
					
					Path baiPath= SamFiles.findIndex(bamPath);
					VirtualFile virtBai = new VirtualFile(baiPath);
					page2content.put("/"+bai_name, virtBai);
					
					final JsonObject idm1 = new JsonObject();
					tracks.add(idm1);
					
					idm1.addProperty("type","AlignmentsTrack");
					idm1.addProperty("trackId",name0+"."+n);
					idm1.addProperty("name",name0);

					final JsonObject idm5 = new JsonObject();
					idm1.add("adapter",idm5);
					idm5.addProperty("type","BamAdapter");

					idm5.add("bamLocation",createUriLocation(bam_name));

					final JsonObject idm10 = new JsonObject();
					idm5.add("index",idm10);

					idm10.add("location",createUriLocation(bai_name));
					idm10.addProperty("indexType","BAI");

					final JsonObject idm15 = new JsonObject();
					idm5.add("sequenceAdapter",idm15);
					idm15.addProperty("type","IndexedFastaAdapter");

					final JsonObject idm17b = new JsonObject();
					idm15.add("fastaLocation",idm17b);
					idm17b.addProperty("uri","rotavirus_rf.fa");
					idm17b.addProperty("locationType","UriLocation");

					final JsonObject idm20b = new JsonObject();
					idm15.add("faiLocation",idm20b);
					idm20b.addProperty("uri","rotavirus_rf.fa.fai");
					idm20b.addProperty("locationType","UriLocation");

					final JsonArray idm23 = new JsonArray();
					idm1.add("assemblyNames",idm23);
					idm23.add(getBuildName());
					}
				return true;
				}
			catch(Exception err) {
				throw err;
				}
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {

		
		try {
			final Map<String,AbstractContent> page2content=new HashMap<String,AbstractContent>() {
				@Override
				public AbstractContent put(String key, AbstractContent value) {
					if(this.containsKey(key)) throw new IllegalArgumentException("duplicate key"+key);
					return super.put(key, value);
					}
				};
			final URL url = new URL(this.jbrowse2url);
			LOG.info(this.jbrowse2url);
			try(InputStream is=url.openStream()) {
				try(ZipInputStream zipin= new ZipInputStream(is)) {
					for(ZipEntry ze=zipin.getNextEntry();ze!=null;ze=zipin.getNextEntry()) {
						if(ze.isDirectory()) continue;
						final String fname=ze.getName();
						if(fname.endsWith(FileExtensions.BAM)) continue;
						if(fname.endsWith(FileExtensions.BAI_INDEX)) continue;
						if(fname.endsWith(FileExtensions.TABIX_INDEX)) continue;
						if(fname.endsWith(FileExtensions.COMPRESSED_VCF)) continue;
						LOG.info(fname);
						try(ByteArrayOutputStream baos = new ByteArrayOutputStream()) {
							IOUtils.copyTo(zipin, baos);
							final byte[] data= baos.toByteArray();
							page2content.put("/"+fname, new FileInZip("/"+fname,data));
							}
						}
					}
				}	
			//add fasta
			final String buildName= IOUtils.getFilenameWithoutCommonSuffixes(this.reference); 
			final String fastaPath="/"+this.reference.getFileName().toString();
			VirtualFile virt = new VirtualFile(this.reference);
			page2content.put(fastaPath, virt);
			
			//add fai
			final Path fai=ReferenceSequenceFileFactory.getFastaIndexFileName(reference);
			final String faiPath="/"+ fai.getFileName().toString();
			virt = new VirtualFile(fai);
			page2content.put(faiPath, virt);
			
			//add dict
			final Path dictf =ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence(this.reference);
			virt = new VirtualFile(dictf);
			page2content.put("/"+ dictf.getFileName().toString(), virt);
			
			
			final JsonObject configJson = new JsonObject();
			final JsonArray assemblies = new JsonArray();
			configJson.add("assemblies", assemblies);
			final JsonObject as=new JsonObject();
			assemblies.add(as);
			as.addProperty("name", buildName);
			final JsonObject sequence = new JsonObject();
			as.add("sequence", sequence);
			sequence.addProperty("type","ReferenceSequenceTrack");
			sequence.addProperty("trackId", buildName+"ReferenceSequenceTrack");
			final JsonObject adapter = new JsonObject();
			sequence.add("adapter",adapter);
			adapter.addProperty("type", "IndexedFastaAdapter");
			JsonObject o2 = new JsonObject();
			o2.addProperty("uri",fastaPath);
			o2.addProperty("locationType","UriLocation");
			adapter.add("fastaLocation",o2);
			o2 = new JsonObject();
			o2.addProperty("uri",faiPath);
			o2.addProperty("locationType","UriLocation");
			adapter.add("faiLocation",o2);
			
			
			
			
			
			final JsonObject idm16 = new JsonObject();
			configJson.add("configuration",idm16);

			final JsonArray idm17 = new JsonArray();
			configJson.add("connections",idm17);

			final JsonObject idm18 = new JsonObject();
			configJson.add("defaultSession",idm18);
			idm18.addProperty("name","My Session");

			final JsonArray tracks = new JsonArray();
			configJson.add("tracks",tracks);
			
			

			for(final String filename: IOUtils.unrollStrings(args)) {
				
				if(filename.endsWith(FileExtensions.BAM)) {
					final SamReaderFactory srf = SamReaderFactory.makeDefault();
					srf.referenceSequence(this.reference);
					try {
						Path bamPath =Paths.get(filename);
						try(SamReader sr=srf.open(bamPath)) {
							if(!sr.hasIndex()) {
								LOG.error("index missing");
								return -1;
								}
							sr.hasIndex();
							String name0 = sr.getFileHeader().getReadGroups().stream().
									map(RG->RG.getSample()).
									filter(S->!StringUtils.isBlank(S)).
									findFirst().
									orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
							String bam_name=name0+FileExtensions.BAM;
							String bai_name=name0+FileExtensions.BAM+FileExtensions.BAI_INDEX;
							int n=1;
							for(;;) {
								if(!page2content.containsKey(bam_name) && !page2content.containsKey( bai_name)) break;
								n++;
								bam_name=name0+"."+n+FileExtensions.BAM;
								bai_name=name0+"."+n+FileExtensions.BAM+FileExtensions.BAI_INDEX;
								}
							VirtualFile virtBam = new VirtualFile(bamPath);
							page2content.put("/"+bam_name, virtBam);
							
							Path baiPath= SamFiles.findIndex(bamPath);
							VirtualFile virtBai = new VirtualFile(baiPath);
							page2content.put("/"+bai_name, virtBai);
							
							final JsonObject idm1 = new JsonObject();
							tracks.add(idm1);
							
							idm1.addProperty("type","AlignmentsTrack");
							idm1.addProperty("trackId",name0+"."+n);
							idm1.addProperty("name",name0);

							final JsonObject idm5 = new JsonObject();
							idm1.add("adapter",idm5);
							idm5.addProperty("type","BamAdapter");

							final JsonObject idm7 = new JsonObject();
							idm5.add("bamLocation",idm7);
							idm7.addProperty("uri",bam_name);
							idm7.addProperty("locationType","UriLocation");

							final JsonObject idm10 = new JsonObject();
							idm5.add("index",idm10);

							final JsonObject idm11 = new JsonObject();
							idm10.add("location",idm11);
							idm11.addProperty("uri",bai_name);
							idm11.addProperty("locationType","UriLocation");
							idm10.addProperty("indexType","BAI");

							final JsonObject idm15 = new JsonObject();
							idm5.add("sequenceAdapter",idm15);
							idm15.addProperty("type","IndexedFastaAdapter");

							final JsonObject idm17b = new JsonObject();
							idm15.add("fastaLocation",idm17b);
							idm17b.addProperty("uri","rotavirus_rf.fa");
							idm17b.addProperty("locationType","UriLocation");

							final JsonObject idm20b = new JsonObject();
							idm15.add("faiLocation",idm20b);
							idm20b.addProperty("uri","rotavirus_rf.fa.fai");
							idm20b.addProperty("locationType","UriLocation");

							final JsonArray idm23 = new JsonArray();
							idm1.add("assemblyNames",idm23);
							idm23.add(buildName);

							}
						}
					catch(Throwable err) {
						LOG.error(err);
						return -1;
						}
					
					}
				else if(filename.endsWith(FileExtensions.COMPRESSED_VCF)) {
					try {
						Path vcfPath =Paths.get(filename);
						Path tbiPath =Paths.get(filename+FileExtensions.COMPRESSED_VCF_INDEX);
						try(VCFFileReader r= new VCFFileReader(vcfPath,true)){
							final VCFHeader header=r.getHeader();
							String name0 = IOUtils.getFilenameWithoutCommonSuffixes(vcfPath);
							String vcf_name=name0+FileExtensions.COMPRESSED_VCF;
							String tbi_name=name0+FileExtensions.COMPRESSED_VCF+FileExtensions.COMPRESSED_VCF_INDEX;
							int n=1;
							for(;;) {
								if(!page2content.containsKey(vcf_name) && !page2content.containsKey( tbi_name)) break;
								n++;
								vcf_name=name0+"."+n+FileExtensions.COMPRESSED_VCF;
								tbi_name=name0+"."+n+FileExtensions.COMPRESSED_VCF+FileExtensions.COMPRESSED_VCF_INDEX;
								}
							VirtualFile virtBam = new VirtualFile(vcfPath);
							page2content.put("/"+vcf_name, virtBam);
							
							VirtualFile virtBai = new VirtualFile(tbiPath);
							page2content.put("/"+tbi_name, virtBai);
							
							final JsonObject idm1 = new JsonObject();
							tracks.add(idm1);
							
							idm1.addProperty("type","VariantTrack");
							idm1.addProperty("trackId",name0+"."+n);
							idm1.addProperty("name",name0);

							final JsonObject idm5 = new JsonObject();
							idm1.add("adapter",idm5);
							idm5.addProperty("type","VcfTabixAdapter");

							final JsonObject idm7 = new JsonObject();
							idm5.add("vcfGzLocation",idm7);
							idm7.addProperty("uri",vcf_name);
							idm7.addProperty("locationType","UriLocation");

							final JsonObject idm10 = new JsonObject();
							idm5.add("index",idm10);

							final JsonObject idm11 = new JsonObject();
							idm10.add("location",idm11);
							idm11.addProperty("uri",tbi_name);
							idm11.addProperty("locationType","UriLocation");
							idm10.addProperty("indexType","TBI");


							final JsonArray idm23 = new JsonArray();
							idm1.add("assemblyNames",idm23);
							idm23.add(buildName);

							}
						}
					catch(Throwable err) {
						LOG.error(err);
						return -1;
						}
					
					
					}
				else
					{
					LOG.warning("skipping "+filename);
					}
				}
			
			page2content.put("/config.json",new JsonConfigFile(configJson));

			final Server server = new Server(this.serverPort);
			final ServletContextHandler context = new ServletContextHandler();
			final JBrowse2Servlet vs = new JBrowse2Servlet(page2content);
			

			
			final ServletHolder sh =new ServletHolder(vs);
	        context.addServlet(sh,"/*");
	        context.setContextPath("/");
	        context.setResourceBase(".");
	        server.setHandler(context);
			
	        final Runnable stop = ()->{
                try {
                    LOG.info("Shutting down ...");

                    sh.stop();
                    server.stop();
                	}
                catch (Exception e) {
                    Thread.currentThread().interrupt();
                    LOG.error(e);
                	}
	        
	        	};
	        
	        Runtime.getRuntime().addShutdownHook(new Thread(stop));
		    
		    LOG.info("Starting server "+getProgramName()+" on port "+this.serverPort);
		    server.start();
		    LOG.info("Server started. Press Ctrl-C to stop. Check your proxy settings ."
		    		+ " Open a web browser at http://localhost:"+this.serverPort+"/jbrowse2 .");
		    server.join();
		    stop.run();
		    return 0;
			}
		catch (final Throwable err) {
			LOG.error(err);
			return -1;
			}
		
		}	

	
	public static void main(String[] args) {
		new JBrowse2Server().instanceMainWithExit(args);
	}

}
