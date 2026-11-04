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

package com.github.lindenb.jvarkit.tools.igvserver;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.locatable.SimpleInterval;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.ContentType;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;


@Program(
		name="igvserver",
		description="Web Server displaying SAM/BAM file. A web interface for jvarkit:tview",
		keywords={"sam","bam","table","visualization","server","web"},
		creationDate = "20261104",
		modificationDate = "20261104"
		)
public class IGVServer extends Launcher {
		private static final Logger LOG=Logger.of(IGVServer.class);
		private static int ID_GENERATOR=0;
		@Parameter(names={"-P","--port","-port"},description="Server listening port")
		private int port=8080;	
		@Parameter(names={"--igv-js-url"},description="Where to download IGV js")
		private String igv_js_url = "https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js";
		@Parameter(names={"--shutdown-after"},description="Stop the server after 'x' seconds.")
		private long shutdownAferSeconds=-1L;


		public static class MyServlet extends HttpServlet {
			private static final long serialVersionUID = 1L;
			private final List<LocalFile> htsFiles;
			private final String igv_js_url;
			MyServlet(
					String igv_js_url,
					List<LocalFile> htsFiles) {
				this.igv_js_url=igv_js_url;
				this.htsFiles = Collections.unmodifiableList(htsFiles);
				}
			@Override
			protected void service(HttpServletRequest req, HttpServletResponse resp)
					throws ServletException, IOException {
				final Set<LocalFile> selected=new HashSet<>();
				
				final String interval_str = StringUtils.ifBlank(req.getParameter("interval"),"");
				
				for(LocalFile lf: this.htsFiles) {
			        	if(!lf.isBamOrCram()) continue;
			        	final String sel=req.getParameter("cb"+lf.id);
			        	System.err.println("cb"+lf.id+"="+sel);
			        	if(!StringUtils.isBlank(sel) && sel.matches("(true|1|yes)")) {
			        		System.err.println("adding to sel "+selected);
			        		selected.add(lf);
			        		}
				 		}
				if(selected.isEmpty()) {
					for(LocalFile lf: this.htsFiles) {
						if(!lf.isBamOrCram()) continue;
						selected.add(lf);
						break;
						}
					}
				final Set<LocalFile> references =new HashSet<>();
				for(LocalFile sel: selected) {
					references.add(sel.reference);
					}
				System.err.println("selected "+selected);
				
				
		         final String charset = StringUtils.ifBlank(req.getCharacterEncoding(), "UTF-8");
		         resp.setContentType("text/html; charset="+charset.toLowerCase());
		         resp.setCharacterEncoding(charset);
		         
		         PrintWriter pw = resp.getWriter();
    	         try {
    	             final XMLStreamWriter w=XMLOutputFactory.newFactory().createXMLStreamWriter(pw);
    	             w.writeStartElement("html");
    	             
    	             w.writeStartElement("head");
    	             
    	             w.writeEmptyElement("meta");
    	             w.writeAttribute("charset", charset);
    	             
    	             w.writeStartElement("title");
    	             w.writeCharacters(interval_str);
    	             w.writeEndElement();//title
		    	     
    	             w.writeStartElement("script");
    	             w.writeAttribute("src",this.igv_js_url);
    	             w.writeCharacters("");
    	             w.writeEndElement();//script
    	             
    	             w.writeEndElement();//read
    	             
    	             
    	             w.writeStartElement("body");
    	             w.writeStartElement("div");
    	             w.writeStartElement("form");
    	             
    	             {
    	             w.writeStartElement("p");
			        	
    	             w.writeStartElement("label");
		        	w.writeAttribute("for", "interval");
		        	w.writeCharacters("Interval:");
		        	w.writeEndElement();
			        	
    	             
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "text");
		        	w.writeAttribute("name", "interval");
		        	w.writeAttribute("id", "interval");
		        	w.writeAttribute("value",interval_str);
		        	
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "submit");
		        	w.writeAttribute("value","Refresh");

			        	
			        	
			        w.writeEndElement();//p
    	             }
    	         	
    	             
    	            for(LocalFile ref: references) {
    	            	w.writeStartElement("div");
    	            	SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(ref.path);
    	            	SimpleInterval region= new IntervalParser(dict).apply(interval_str).orElse(null);
    	            	if(region!=null) {
    	            		w.writeStartElement("div");
    	            		w.writeAttribute("id", "ref"+ref.id);
    	            		w.writeEndElement();
    	            		
    	            	
	    	            	w.writeStartElement("script");
	    	            	w.writeCharacters("const options_"+ref.id+"={");
	    	            	w.writeCharacters("reference:{");
	    	            	w.writeCharacters("id:\"ref"+ref.id+"\",");
	    	            	w.writeCharacters("name:\""+SequenceDictionaryUtils.getBuildName(dict).orElse("genome") +"\",");
	    	            	w.writeCharacters("fastaURL:\"/files"+ref.path.toAbsolutePath().toString() +"\",");
	    	            	w.writeCharacters("indexURL:\"/files"+ref.index.path.toAbsolutePath().toString() +"\"");
	    	            	w.writeCharacters("}");
	    	            	w.writeCharacters(",locus:\""+region.toNiceString()+"\"");
	
	    	            	
	    	            	w.writeCharacters(",tracks:[");
	    	            	w.writeCharacters(selected.stream().
	    	            		filter(SEL->SEL.reference.equals(ref)).
	    	            		map(S->S.toTrack()).collect(Collectors.joining(","))
	    	            		);
	    	            	
	    	            	w.writeCharacters("]};\n");
	    	            	w.writeCharacters("var igvDiv"+ref.id+" = document.getElementById(\"ref"+ref.id+"\");\n");
	    	            	w.writeCharacters("  igv.createBrowser(igvDiv"+ref.id+", options_"+ref.id+")\n"
	    	            			+ "        .then(function (browser) {\n"
	    	            			+ "            console.log(\"Created IGV browser\")\n"
	    	            			+ "        });\n");
	    	            	
	    	            	
	    	            	w.writeEndElement();//script
	    	            	}
    	            	w.writeEndElement();
    	            	}
    	             
    	            
    	            w.writeEmptyElement("hr");
    	             
    	         
    	             for(LocalFile lf: this.htsFiles) {
    			        	if(!lf.isBamOrCram()) continue;
    			        	w.writeStartElement("p");
    			        	
    			        	w.writeEmptyElement("input");
    			        	w.writeAttribute("type", "checkbox");
    			        	if(selected.contains(lf)) w.writeAttribute("checked", "true");
    			        	w.writeAttribute("name", "cb"+lf.id);
    			        	w.writeAttribute("id", "cb"+lf.id);
    			        	w.writeAttribute("value", "1");
    			        	
    			        	w.writeStartElement("label");
    			        	w.writeAttribute("for", "cb"+lf.id);
    			        	w.writeCharacters(lf.path.getFileName().toString());
    			        	w.writeEndElement();
    			        	
    			        	
    			        	w.writeEndElement();//p
    			        	}
    	             
    	             w.writeEndElement();//form
    	             w.writeEndElement();//div
    	             w.writeEndElement();//body
    	             w.writeEndElement();//html
    	             w.flush();
    	             w.close();
    	         	}
    	         catch(XMLStreamException err) {
    	        	 err.printStackTrace();
    	        	 throw new IOException(err);
    	         	}
		        pw.flush();
		        pw.close();
				}
			
			}
		private static class LocalFile {
			final int id = (++ID_GENERATOR);
			final Path path;
			LocalFile index=null;
			LocalFile reference=null;
			LocalFile(final Path path) {
				IOUtil.assertFileIsReadable(path);
				this.path = path;
				}
			boolean hasSuffix(final String sfx) {
				return path.getFileName().toString().endsWith(sfx);
				}
			boolean isBamOrCram() {
				return hasSuffix(FileExtensions.CRAM) || hasSuffix(FileExtensions.BAM);
				}
			boolean isFasta() {
				return FileExtensions.FASTA.stream().anyMatch(SFX->hasSuffix(SFX));
				}
			public String getURI() {
				return this.path.toAbsolutePath().toString();
				/*
				return  new StringBuilder(RAW_PATH)
						.append("/")
						.append(StringUtils.md5(path.toString()))
						.append("/")
						.append(this.path.getFileName().toString())
						.toString();*/
				}
			@Override
			public int hashCode() {
				return path.hashCode();
				}
			String toTrack() {
				return new StringBuilder("{")
					.append("type:\"alignment\",format:\"")
					.append(hasSuffix(FileExtensions.BAM)?"bam":"cram")
					.append("\",url:\"")
					.append("/files"+this.path.toAbsolutePath().toString())
					.append("\",indexURL:\"")
					.append("/files"+this.index.path.toAbsolutePath().toString())
					.append("\",name:\"")
					.append(IOUtils.getFilenameWithoutCommonSuffixes(this.path))
					.append("\",height:600")
					.append("}")
					.toString();
				}
			
			@Override
			public boolean equals(Object obj) {
				if(obj==this) return true;
				if(obj==null ||!(obj instanceof LocalFile)) return false;
				return this.path.equals(LocalFile.class.cast(obj).path);
 				}
			@Override
			public String toString() {
				return getURI()+" = " +this.path;
				}
			}
		
		
		
		@Override
		public int doWork(final List<String> args) {
			Server server = null;
			try
				{
				final Map<SAMSequenceDictionary,LocalFile> dict2ref=new HashMap<>();
				final List<LocalFile> hts_files= new ArrayList<>();
				for(Path htsFile: IOUtils.unrollPaths(args)) {
					final String fname = htsFile.getFileName().toString();
					if(fname.endsWith(FileExtensions.BAM) || fname.endsWith(FileExtensions.CRAM)) {
						final Path index = SamFiles.findIndex(htsFile);
						if(index==null) {
							LOG.error("cannot find index for "+htsFile);
							return -1;
							}
						final LocalFile lf=new LocalFile(htsFile);
						lf.index=new LocalFile(index);
						hts_files.add(lf);
						hts_files.add(lf.index);
						}
					else if(FileExtensions.FASTA.stream().anyMatch(SFX->fname.endsWith(SFX))) {
						LocalFile lf=new LocalFile(htsFile);
						
						final Path fai = ReferenceSequenceFileFactory.getFastaIndexFileName(lf.path);
						IOUtil.assertFileIsReadable(fai);
						lf.index=new LocalFile(fai);
						
						final SAMSequenceDictionary dict1 = SequenceDictionaryUtils.extractRequired(lf.path);
						
						dict2ref.put(dict1,lf);
						}
					else
						{
						LOG.error("undefined file type "+htsFile);
						return -1;
						}
					}
				
				for(LocalFile lf1: hts_files) {
					if(!lf1.isBamOrCram()) continue;
					final SAMSequenceDictionary dict1 = SequenceDictionaryUtils.extractRequired(lf1.path);
					lf1.reference=dict2ref.entrySet().stream()
							.filter(D->SequenceUtil.areSequenceDictionariesEqual(dict1, D.getKey()))
							.map(D->D.getValue())
							.findAny().orElse(null);
					if(lf1.reference==null) {
						LOG.error("cannot link ref for bam file "+lf1.path+" "+dict2ref+" "+dict1);
						return -1;
						}
					}
				
				
				if(hts_files.isEmpty()) {
					LOG.error("no file was found");
					return -1;
					}
				System.err.println(hts_files);
				
				server = new Server(this.port);
			
		        
		        final ServletContextHandler ctx = new ServletContextHandler();
		        ctx.setContextPath("/");
		        ctx.setResourceBase("/");
				ctx.addServlet(new ServletHolder("my", new MyServlet(
						this.igv_js_url,
						hts_files
						)),"/igv/*");
				
				
		        final ServletHolder files = new ServletHolder("static_file_server", DefaultServlet.class);
		        files.setInitParameter("baseResource","/");
		        files.setInitParameter("resourceBase","/");
		        files.setInitParameter("dirAllowed", "true");   // directory listing; set to false to disable
		        files.setInitParameter("welcomeServlets", "false");
		        files.setInitParameter("pathInfoOnly", "true");
		        ctx.addServlet(files, "/files/*");

		        
		        
		        
		        
		        // 3. Mount it on the server
		        server.setHandler(ctx);
				
		        
				
				
				
		        LOG.info("Starting IGVServer on http://localhost:"+this.port);
				server.start();
				if(this.shutdownAferSeconds>0)
					{
					final Server theServer = server;
					new java.util.Timer().schedule( 
					        new java.util.TimerTask() {
					            @Override
					            public void run() {
					                LOG.info("automatic shutdown after "+shutdownAferSeconds);
					                try {
					                	theServer.stop();
					                	}
					                catch(final Throwable err2) {
					                	LOG.error(err2);
					                	}
					            }
					        }, 
					        1000 * this.shutdownAferSeconds 
							);
					}
				
				server.join();
				return 0;
				}
			catch(final Throwable err)  {
				LOG.error(err);
				return -1;
				}
			finally
				{
				if(server!=null)
					{
					server.destroy();
					}
				}
			}

		public static void main(final String args[]) {
			new IGVServer().instanceMainWithExit(args);
			}
	}
