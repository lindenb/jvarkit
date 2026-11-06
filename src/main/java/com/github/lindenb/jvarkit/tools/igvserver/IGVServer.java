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
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;

import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.BooleanUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.HasName;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.net.ContentType;
import com.github.lindenb.jvarkit.net.UrlSupplier;
import com.github.lindenb.jvarkit.net.UrlSupplier.LabelledUrl;
import com.github.lindenb.jvarkit.samtools.util.IntervalParser;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;

/**
BEGIN_DOC

## Motivation

creates a IGV.js browser that is able to visualize multiple BAM/vcf files


## Input

input is a set of bam/vcf/cram/fasta file. A fasta file compatible with the bam/vcf is required.


## Example:

```
 $ java -jar dist/jvarkit.jar igvserver --bed roi.bed src/test/resources/rotavirus_rf.fa src/test/resources/S*.bam src/test/resources/rotavirus_rf.vcf.gz
 $ java -jar dist/jvarkit.jar igvserver igv_sessions1.xml igv_sessions2.xml 
 
$ find /path/to -type  f -name "*.cram" -o -name "*.fa" > input.list
$ java -jar dist/jvarkit.jar igvserver input.list 
 
```


END_DOC
 */
@Program(
		name="igvserver",
		description="Web Server displaying SAM/BAM/VCF file and using the IGV.js API.",
		keywords={"sam","bam","table","visualization","server","web"},
		creationDate = "20261104",
		modificationDate = "20261106",
		jvarkit_amalgamion = true,
		generate_doc = true
		)
public class IGVServer extends Launcher {
		private static final Logger LOG=Logger.of(IGVServer.class);
		private static final String STATIC_PATH="static";
		private static final String PARAM_INTERVAL="interval";
		private static final String PARAM_TRACK_HEIGHT="th";
		private static final int MIN_TRACK_HEIGHT = 100;
		private static boolean FIRST_CALL=true;
		private static int ID_GENERATOR=0;
		@Parameter(names={"-P","--port","-port"},description="Server listening port")
		private int port=8080;	
		@Parameter(names={"--igv-js-url"},description="Where to download IGV js")
		private String igv_js_url = "https://cdn.jsdelivr.net/npm/igv@3.8.10/dist/igv.min.js";
		@Parameter(names={"--shutdown-after"},description="Stop the server after 'x' seconds.")
		private long shutdownAferSeconds=-1L;
		@Parameter(names={"--bed","--roi"},description="Optional Bed of regions of interest displated as a html 'select' item")
		private Path roiBedPath=null;

		
		private static class ConfigParam {
			final String id;
			final String name;
			final String label;
			ConfigParam(final String id,final String name, final String label) {
				this.id=id;
				this.label=label;
				this.name=name;
				}
			}

		
		private static class ConfigFlag extends ConfigParam {
			final boolean def;
			ConfigFlag(final String id,final String name, final String label,boolean def) {
				super(id,name,label);
				this.def = def;
				}
			}

		private static final List<ConfigFlag> CONFIG_FLAGS= Arrays.asList(
				new ConfigFlag("showcov"," showCoverage"," Show coverage depth track. ",true),
				new ConfigFlag("showaln"," showAlignments"," Show individual alignments. ",true),
				new ConfigFlag("aspair"," viewAsPairs"," If true, paired reads are drawn connected with a line. ",false),
				new ConfigFlag("pairsupp"," pairsSupported"," If false, mate information in paired reads is ignored during downsampling and the 'View as Pairs' option is removed from the alignment track menu.",false),
				new ConfigFlag("softclip"," showSoftClips"," Show soft-clipped regions ",true),
				new ConfigFlag("mismatches"," showMismatches"," Highlight alignment bases which do not match the reference. ",true),
				new ConfigFlag("allbases"," showAllBases"," Show all bases of the read sequence. ",false),
				new ConfigFlag("instxt"," showInsertionText"," Show number of bases for insertions inline when zoomed in. ",true),
				new ConfigFlag("deltxt"," showDeletionText"," Show number of bases deleted inline when zoomed in. ",true),
				new ConfigFlag("qcfail"," vendorFailed"," filter alignments marked as failing vendor quality checks (bit 0x200) ",true),
				new ConfigFlag("dups"," duplicates"," filter alignments marked as a duplicate (bit 0x400) ",true),
				new ConfigFlag("second"," secondary"," filter alignments marked secondary (bit 0x100) ",false),
				new ConfigFlag("suppl"," supplementary"," filter alignments marked as supplmentary (bit 0x800) ",false)
				);

		
		private static class ConfigEnum extends ConfigParam {
			private final List<String> values;
			ConfigEnum(final String id,final String name, final String label,String values) {
				super(id,name,label);
				this.values = Arrays.asList(CharSplitter.COMMA.split(values));
				}
			}
		
		private static final List<ConfigEnum> CONFIG_ENUMS= Arrays.asList(
				new ConfigEnum("co","colorBy","Color alignment by property","none,strand,firstOfPairStrand,pairOrientation,tlen,unexpectedPair,tag:YC"),
				new ConfigEnum("gb","groupBy","Group alignments by property","none,strand,firstOfPairStrand,pairOrientation,mateChr,chimeric,supplementary,readOrder"),
				new ConfigEnum("dm","displayMode","Display Mode","EXPANDED,SQUISHED,FULL")
			);

		
		public static class IGVServlet extends HttpServlet {
			private static final long serialVersionUID = 1L;
			private final List<LocalFile> htsFiles;
			private final String igv_js_url;
			private final List<Interval> roi_interval;

			IGVServlet(
					String igv_js_url,
					List<LocalFile> htsFiles,
					List<Interval> roi_interval
					) {
				this.igv_js_url=igv_js_url;
				this.htsFiles = Collections.unmodifiableList(htsFiles);
				this.roi_interval = Collections.unmodifiableList(roi_interval);
				}
			
			private final boolean parseBoolean(String s) {
				return BooleanUtils.parseBoolean(s).orElse(Boolean.FALSE);
				}
			
			private String resolveStatic(final Path path) {
				return STATIC_PATH+ path.toAbsolutePath().toString();
				}
			
			private String toBamTrack(LocalFile bam,int height,
					final Map<String,Boolean> configFlags,
					final Map<String,String> configEnums
					) {
				final StringBuilder sb= new StringBuilder("{");
				if(bam.isBamOrCram()) {
					sb
						.append("type:\"alignment\",format:\"")
						.append(bam.hasSuffix(FileExtensions.BAM)?"bam":"cram")
						.append("\",url:\"")
						.append("/"+resolveStatic(bam.path))
						.append("\",indexURL:\"")
						.append("/"+resolveStatic(bam.index.path))
						.append("\",name:\"")
						.append(bam.getName())
						.append("\",height:")
						.append(Math.max(height, MIN_TRACK_HEIGHT))
						;
					}
				else if(bam.isVcf()) {
					sb
						.append("type:\"variant\",format:\"vcf\",")
						.append("url:\"")
						.append("/"+resolveStatic(bam.path))
						.append("\",indexURL:\"")
						.append("/"+resolveStatic(bam.index.path))
						.append("\",name:\"")
						.append(bam.getName())
						.append("\"")
						;
					}
				for(final String k: configFlags.keySet()) {
					sb.append(",").append(k).append(":").append(configFlags.get(k));
					}
			
				for(final String k: configEnums.keySet()) {
					sb.append(",").append(k).append(":\"").append(configEnums.get(k)).append("\"");
					}
				
				return sb.append("}").toString();
				}
			
			private Locatable  parseInterval(final SAMSequenceDictionary dict, String interval_str) {
				if(StringUtils.isBlank(interval_str)) return null;
				interval_str=interval_str.replaceAll(",", "");
				if(StringUtils.isBlank(interval_str)) return null;
				String[] tokens= interval_str.split("\\s+");
				if(tokens.length==3 && StringUtils.isInteger(tokens[1]) && StringUtils.isInteger(tokens[2])) {
					parseInterval(dict,tokens[0]+":"+tokens[1]+"-"+tokens[2]);
					}
				if(tokens.length==2 && StringUtils.isInteger(tokens[1])) {
					parseInterval(dict,tokens[0]+":"+tokens[1]+"-"+tokens[1]);
					}
						
				return new IntervalParser(dict).apply(interval_str).orElse(null);
				}
			
			private String toString(Locatable l) {
				return l.getContig()+":"+l.getStart()+"-"+l.getEnd();
				}
			private String toNiceString(Locatable l) {
				return l.getContig()+":"+StringUtils.niceInt(l.getStart())+"-"+StringUtils.niceInt(l.getEnd());
				}
			@Override
			protected void service(final HttpServletRequest req, HttpServletResponse resp)
					throws ServletException, IOException {
				final List<LocalFile> bam_vcf_files = this.htsFiles
						.stream()
						.filter(F->F.isBamOrCramOrVcf())
						.collect(Collectors.toList());
				final int n_random = 10;//number of random files we can choose
				final Set<LocalFile> selected=new HashSet<>();
				final Map<String,Boolean> configFlags = new HashMap<String, Boolean>(CONFIG_FLAGS.size());
				final Map<String,String> configEnums = new HashMap<String,String>(CONFIG_ENUMS.size());
				final String interval_str = StringUtils.ifBlank(req.getParameter(PARAM_INTERVAL),
						this.roi_interval.isEmpty()?
						"":
						toString(this.roi_interval.get(0))
						);
				for(ConfigFlag flag: CONFIG_FLAGS) {
					configFlags.put(flag.name, FIRST_CALL?flag.def: parseBoolean(req.getParameter(flag.id)));
					}
				for(ConfigEnum cfg: CONFIG_ENUMS) {
					configEnums.put(cfg.name, FIRST_CALL?cfg.values.get(0): StringUtils.ifBlank(req.getParameter(cfg.id),cfg.values.get(0)));
					}
				
				final int track_height = (FIRST_CALL?300:Math.max(StringUtils.parseInt(req.getParameter(PARAM_TRACK_HEIGHT)).orElse(0),MIN_TRACK_HEIGHT));
				
				for(LocalFile lf: bam_vcf_files) {
			        	final String sel=req.getParameter(lf.getAttribute());
			        	if(!StringUtils.isBlank(sel) && sel.matches("(true|1|yes)")) {
			        		selected.add(lf);
			        		}
				 		}
				if(selected.isEmpty() && FIRST_CALL) {
					for(LocalFile lf: bam_vcf_files) {
						selected.add(lf);
						break;
						}
					}
				final Set<LocalFile> references =new HashSet<>();
				for(LocalFile sel: selected) {
					references.add(sel.reference);
					}
				
				
		         final String charset = StringUtils.ifBlank(req.getCharacterEncoding(), "UTF-8");
		         resp.setContentType( ContentType.TEXT_XHTML.getMimeType()+"; charset="+charset.toLowerCase());
		         resp.setCharacterEncoding(charset);
		         
		         final PrintWriter pw = resp.getWriter();
    	         try {
    	             final XMLStreamWriter w=XMLOutputFactory.newFactory().createXMLStreamWriter(pw);
    	             //w.writeStartDocument(charset, "1.0");
    	             w.writeDTD("<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">");
    	             w.writeStartElement("html");
    	             w.writeAttribute("xml:lang", "en");
    	             w.writeAttribute("lang", "en");
    	             w.writeAttribute("xmlns", "http://www.w3.org/1999/xhtml");
    	             
    	             
    	             w.writeStartElement("head");
    	             
    	             
    	             w.writeEmptyElement("meta");
    	             w.writeAttribute("charset", charset);
    	             w.writeEmptyElement("meta");
    	             w.writeAttribute("http-equiv", "author");
    	             w.writeAttribute("content", "Pierre Lindenbaum");

    	             
    	             w.writeStartElement("title");
    	             w.writeCharacters(interval_str);
    	             w.writeEndElement();//title
		    	     
    	             w.writeStartElement("script");
    	             w.writeAttribute("src",this.igv_js_url);
    	             w.writeEndElement();//script
    	             
    	             
    	             w.writeCharacters("\n");
    	             w.writeStartElement("script");
    	             w.writeAttribute( "type","text/javascript");
    	             // list of all sample IDs
    	             w.writeCharacters("var hts_ids=["+bam_vcf_files.stream().map(F->"\""+F.getAttribute()+"\"").collect(Collectors.joining(","))+"];\n");
    	             // function to toggle the samples
    	             w.writeCharacters(
    	            		"function toggle_checkboxes(a,b) {\n"
     	             		+ "for (let i = 0; i < a.length; i++)  {document.getElementById(a[i]).checked=b;}\n"
     	             		+ "}\n");
    	             // select sample by name and choose n_random samples
    	             w.writeCharacters("function sel_random(name) {\n"
    	             		+ " toggle_checkboxes(hts_ids,false);\nvar a=[];\nif(name!=null) a.push(name);\n"
    	            		+ " for(var i=0; i< "+n_random+" && i< hts_ids.length ;i++) { a.push( hts_ids[ Math.floor( Math.random()*hts_ids.length ) ] ); }\n"
    	            		+ " toggle_checkboxes(a,true);\n"
    	             		+ "}\n");
    	             
    	             w.writeCharacters("\nwindow.addEventListener('load', function() {\n"
      	            		+ "  var sel =  document.getElementById(\"roi\"); \n"
     	            		+ "  var input1 =  document.getElementById(\""+PARAM_INTERVAL+"\"); \n"
     	            		+ "  var form1 =  document.getElementById(\"form\"); \n"
    	            		+ "  sel.addEventListener(\"change\", function(event2)  {if(sel.value){ input1.value=sel.value; form1.submit();} });\n"
    	             		+ "});\n");
    	             
    	             
    	             w.writeEndElement();//script
    	             
    	             
    	             w.writeStartElement("style");

                     w.writeCharacters(
                             "body {background-color:#f0f0f0;color:#070707;font-size:18px;}"
                             + "h1 {text-align:center;color:#070707;}"
                             + ".grid-container {display: grid; grid-template-columns: auto auto auto auto auto auto;grid-gap: 10px;  padding: 10px;}"
                             + ".comment {background-color:khaki;text-align:center;font-size:14px;}"
                             + ".highlight {background-color:#DB7093;}"
                             + ".message {color:red;}"
                             + "span.small {font-size: smaller;}"
                             );
                     w.writeEndElement();//style
    	             
    	             w.writeEndElement();//head
    	             
    	             
    	             w.writeStartElement("body");
    	             w.writeStartElement("div");
    	             
    	             w.writeStartElement("h1");
    	             w.writeCharacters(interval_str);
    	             w.writeEndElement();//h1
    	             
    	             w.writeStartElement("form");
			         w.writeAttribute("id", "form");
			         //too many files to handle GET, use POST
			         if(bam_vcf_files.size()>100) {
			        	w.writeAttribute("method", "POST"); 
			         	}
    	             
    	             {
    	             w.writeStartElement("p");
    	             
    	             if(!this.roi_interval.isEmpty()) {
    	            	w.writeStartElement("label");
			        	w.writeAttribute("for", "roi");
			        	w.writeCharacters("ROI:");
			        	w.writeEndElement();
    	            	 
    	            	w.writeStartElement("select");
    	            	w.writeAttribute("id", "roi");
			        	for(Interval roi: this.roi_interval) {
			        		w.writeStartElement("option");
			        		if(toString(roi).equals(interval_str)) {
			        			w.writeAttribute("selected", "true");
			        			}
			        		w.writeAttribute("value", toString(roi));
			        		w.writeCharacters(toNiceString(roi));
			        		w.writeEndElement();
			        		}
    	            	w.writeEndElement();
    	            	
			        	w.writeCharacters(" ");
    	             } 
    	             
    	             
			        // INTERVAL
    	            w.writeStartElement("label");
		        	w.writeAttribute("for", PARAM_INTERVAL);
		        	w.writeCharacters("Interval:");
		        	w.writeEndElement();
			        	
    	             
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "text");
		        	w.writeAttribute("name", PARAM_INTERVAL);
		        	w.writeAttribute("id", PARAM_INTERVAL);
		        	w.writeAttribute("value",interval_str);
		        	
		        	
		        	w.writeCharacters(" ");
		        	
		        	// TRACK HEIGHT
		        	
		        	w.writeStartElement("label");
		        	w.writeAttribute("for", PARAM_TRACK_HEIGHT);
		        	w.writeCharacters("Track Height:");
		        	w.writeEndElement();	
    	             
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "text");
		        	w.writeAttribute("name", PARAM_TRACK_HEIGHT);
		        	w.writeAttribute("id", PARAM_TRACK_HEIGHT);
		        	w.writeAttribute("value",String.valueOf(track_height));
		        	
		        	w.writeCharacters(" ");
		        	
		        	for(ConfigFlag flag: CONFIG_FLAGS) {
			        	w.writeStartElement("span");
			        	w.writeAttribute("class", "small");

			        	w.writeEmptyElement("input");
			        	w.writeAttribute("type", "checkbox");
			        	if((FIRST_CALL && flag.def) || parseBoolean(req.getParameter(flag.id))) w.writeAttribute("checked", "true");
			        	w.writeAttribute("name", flag.id);
			        	w.writeAttribute("id", flag.id);
			        	w.writeAttribute("value", "1");
			        	
			        	w.writeStartElement("label");
			        	w.writeAttribute("for",flag.id);
			        	w.writeCharacters(flag.label);
			        	w.writeEndElement();
			        	
			        	w.writeEndElement();//span
			        	w.writeCharacters(" ");
			        	}
		        	
		        	for(ConfigEnum flag: CONFIG_ENUMS) {
			        	w.writeStartElement("span");
			        	w.writeAttribute("class", "small");

			        	
			        	w.writeStartElement("label");
			        	w.writeAttribute("for",flag.id);
			        	w.writeCharacters(flag.label);
			        	w.writeEndElement();
			        	
			        	
			        	w.writeStartElement("select");
			        	w.writeAttribute("name", flag.id);
			        	w.writeAttribute("id", flag.id);
			        	for(int x=0;x< flag.values.size();++x) {
			        		w.writeStartElement("option");
				        	if((x==0 && FIRST_CALL) || flag.values.get(x).equals(req.getParameter(flag.id))) w.writeAttribute("selected", "true");
				        	w.writeCharacters(flag.values.get(x));
				        	w.writeEndElement();//option
			        		}
			        	w.writeEndElement();//select
			        	
			        	w.writeEndElement();//span
			        	w.writeCharacters(" ");
			        	}

		        	
		        	
		        	
		        	w.writeCharacters(" ");
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "submit");
		        	w.writeAttribute("value","Refresh");

			        	
			        	
			        w.writeEndElement();//p
    	             }
    	         	
    	             
    	            for(LocalFile ref: references) {
    	            	w.writeStartElement("div");
    	            	final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(ref.path);
    	            	final Locatable region= parseInterval(dict,interval_str);//
    	            	if(region==null) {
    	            		if(!StringUtils.isBlank(interval_str)) {
	    	            		w.writeStartElement("p");
	    	            		w.writeAttribute("class", "message");
	    	            		w.writeCharacters("cannot parse interval "+interval_str);
	    	            		w.writeEndElement();
	    	            		}
    	            		}
    	            	else {
    	            		w.writeStartElement("div");
    	            		w.writeAttribute("id", "ref"+ref.id);
    	            		w.writeEndElement();
    	            		
    	            	
    	            		
    	            		w.writeStartElement("div");
    	            		w.writeStartElement("table");
    	            		w.writeStartElement("tr");
    	            		w.writeStartElement("th");
    	            		w.writeCharacters("Database");
    	            		w.writeEndElement();//th
    	            		w.writeStartElement("th");
    	            		w.writeCharacters("URL");
    	            		w.writeEndElement();//th

    	            		
    	            		w.writeEndElement();//tr
    	            		for(LabelledUrl url: new UrlSupplier(dict).of(region)) {
    	            			w.writeStartElement("tr");
    	            			
    	            			w.writeStartElement("td");
        	            		w.writeCharacters(url.getLabel());
        	            		w.writeEndElement();//td
        	            		
    	            			w.writeStartElement("td");
    	            			w.writeStartElement("a");
        	            		w.writeAttribute("title", url.getLabel());
        	            		w.writeAttribute("href", url.getUrl());
        	            		w.writeAttribute("target", "_blank");
        	            		w.writeCharacters(url.getUrl());
        	            		w.writeEndElement();//a
        	            		w.writeEndElement();//td

        	            		
    	            			w.writeEndElement();//tr
    	            			}
    	            		w.writeEndElement();//table
    	            		w.writeEndElement();//div
    	            		
	    	            	w.writeStartElement("script");
	    	            	w.writeCharacters("const options_"+ref.id+"={");
	    	            	w.writeCharacters("reference:{");
	    	            	w.writeCharacters("id:\"ref"+ref.id+"\",");
	    	            	w.writeCharacters("name:\""+SequenceDictionaryUtils.getBuildName(dict).orElse("genome") +"\",");
	    	            	w.writeCharacters("fastaURL:\"/"+resolveStatic(ref.path) +"\",");
	    	            	w.writeCharacters("indexURL:\"/"+resolveStatic(ref.index.path) +"\"");
	    	            	w.writeCharacters("}");
	    	            	w.writeCharacters(",locus:\""+toString(region)+"\"");
	
	    	            	
	    	            	w.writeCharacters(",tracks:[");
	    	            	w.writeCharacters(selected.stream().
	    	            		filter(SEL->SEL.reference.equals(ref))
	    	            		.map(S->toBamTrack(S,track_height,configFlags,configEnums))
	    	            		.collect(Collectors.joining(","))
	    	            		);
	    	            	
	    	            	w.writeCharacters("]};\n");
	    	            	w.writeCharacters("var igvDiv"+ref.id+" = document.getElementById(\"ref"+ref.id+"\");\n"
	    	            			+ "  igv.createBrowser(igvDiv"+ref.id+", options_"+ref.id+")\n"
	    	            			+ "        .then(function (browser) {\n"
	    	            			+ "            console.log(\"Created IGV browser\");\n"
	    	            			+ "             browser.on('locuschange', function (referenceFrameIdx) {\n"
	    	            			+ "                     if(browser.referenceFrameList==undefined || browser.referenceFrameList.length==0) return;\n"
	    	            			+ "                     var referenceFrame= browser.referenceFrameList[0];"
	    	            			+ "                     document.getElementById(\""+PARAM_INTERVAL+"\").value=referenceFrame.chr + \":\" + parseInt(referenceFrame.start) + \"-\" + parseInt(referenceFrame.end);\n"
	    	            			+ "                     });\n"
	    	            			+ "        });\n");
	    	            	
	    	            	
	    	            	w.writeEndElement();//script
	    	            	}
    	            	
    	            	w.writeEndElement();
    	            	}
   	            

    	            
    	            w.writeEmptyElement("hr");
    	             
    	            // button to select ALL
		        	w.writeStartElement("button");
		        	w.writeAttribute("type", "button");
		        	w.writeAttribute("title", "Select all Samples");
		        	w.writeAttribute("onclick", "toggle_checkboxes(hts_ids,true)");
		        	w.writeCharacters("Select All");
		        	w.writeEndElement();
		        	w.writeCharacters(" ");
    	            
    	            // button to select None
		        	w.writeStartElement("button");
		        	w.writeAttribute("type", "button");
		        	w.writeAttribute("title", "Select no Samples");
		        	w.writeAttribute("onclick", "toggle_checkboxes(hts_ids,false)");
		        	w.writeCharacters("Select None");
		        	w.writeEndElement();
		        	w.writeCharacters(" ");
		        	
    	            // button to select Random
		        	w.writeStartElement("button");
		        	w.writeAttribute("type", "button");
		        	w.writeAttribute("title", "Select "+n_random+" Samples");
		        	w.writeAttribute("onclick", "sel_random(null)");
		        	w.writeEntityRef("#x1F3B2");
		        	w.writeCharacters(" ");
		        	w.writeCharacters("Random "+n_random);
		        	w.writeEndElement();
		        	w.writeCharacters(" ");
 	
		        	 // button to submit
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "submit");
		        	w.writeAttribute("value","Refresh");
    	            
    	            w.writeStartElement("div");
    	            w.writeAttribute("class", "grid-container");
    	            for(LocalFile lf: bam_vcf_files) {
    			        	w.writeStartElement("div");
    			        	
    			        	w.writeEmptyElement("input");
    			        	w.writeAttribute("type", "checkbox");
    			        	if(selected.contains(lf)) w.writeAttribute("checked", "true");
    			        	w.writeAttribute("name", lf.getAttribute());
    			        	w.writeAttribute("id", lf.getAttribute());
    			        	w.writeAttribute("value", "1");
    			        	
    			        	w.writeStartElement("label");
    			        	w.writeAttribute("for",lf.getAttribute());
    			        	w.writeAttribute("title",lf.path.toString());
    			        	w.writeCharacters(lf.getName());
    			        	w.writeEndElement();
    			        	
    			        	if(bam_vcf_files.size()  > n_random) {
	    			        	w.writeCharacters(" ");
	    			        	w.writeStartElement("button");
	    			        	w.writeAttribute("type", "button");
	    			        	w.writeAttribute("title", "choose "+n_random+" files");
	    			        	w.writeAttribute("onclick", "sel_random(\""+lf.getAttribute()+"\")");
	    			        	w.writeEntityRef("#x1F3B2");
	    			        	w.writeCharacters(" ");
	    			        	w.writeCharacters(String.valueOf(n_random));
	    			        	w.writeEndElement();
	    			        	}
    			        	
    			        	w.writeEndElement();//div
    			        	}
    	             w.writeEndElement();
    	             
    	             
		        	w.writeEmptyElement("p");
		        	w.writeEmptyElement("input");
		        	w.writeAttribute("type", "submit");
		        	w.writeAttribute("value","Refresh");
    	             
    	             
    	             w.writeEndElement();//form
    	             w.writeEndElement();//div
    	             
                     w.writeEmptyElement("hr");
                     w.writeStartElement("div");
                     w.writeCharacters("Author: Pierre Lindenbaum. ");
                     w.writeCharacters(JVarkitVersion.getInstance().getLabel());
                     w.writeEndElement();

    	             
    	             w.writeEndElement();//body
    	             w.writeEndElement();//html
    	             //w.writeEndDocument();
    	             w.flush();
    	             w.close();
    	         	}
    	         catch(XMLStreamException err) {
    	        	 err.printStackTrace();
    	        	 throw new IOException(err);
    	         	}
		        pw.flush();
		        pw.close();
				FIRST_CALL=false;
				}
			}
		private static class LocalFile implements HasName{
			final int id = (++ID_GENERATOR);
			final Path path;
			/** optional index */
			LocalFile index=null;
			/** optional reference */
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
			boolean isVcf() {
				return hasSuffix(FileExtensions.COMPRESSED_VCF) || hasSuffix(FileExtensions.COMPRESSED_VCF_BGZ);
				}
			boolean isBamOrCramOrVcf() {
				return isVcf() || isBamOrCram();
				}
			
			
			@Override
			public int hashCode() {
				return path.hashCode();
				}
			
			public String getAttribute() {
				return "cb"+this.id;
			}
			
			@Override
			public String getName() {
				return IOUtils.getFilenameWithoutCommonSuffixes(this.path);
				}
			@Override
			public boolean equals(Object obj) {
				if(obj==this) return true;
				if(obj==null ||!(obj instanceof LocalFile)) return false;
				return this.path.equals(LocalFile.class.cast(obj).path);
 				}
			@Override
			public String toString() {
				return String.valueOf(this.path);
				}
			}
		
		
		
		@Override
		public int doWork(final List<String> args) {
			Server server = null;
			try
				{
				final Map<SAMSequenceDictionary,LocalFile> dict2ref=new HashMap<>();
				final List<LocalFile> hts_files= new ArrayList<>();
				final List<Path> input_files = IOUtils.unrollPaths(args);
				final List<Interval> roi_intervals = new Vector<>();
				
				// input is a set of igv sessions XML file
				if(input_files.stream().allMatch(F->F.getFileName().toString().endsWith(".xml"))) {
					final List<Path> input_files2= new ArrayList<Path>();
					final DocumentBuilder db= DocumentBuilderFactory.newInstance().newDocumentBuilder();
					for(Path xmlPath: input_files) {
						// parse igv sessions
						final Document dom= db.parse(xmlPath.toFile());
						final Element root = dom.getDocumentElement();
						if(root==null || !"Session".equals(root.getNodeName())) {
							LOG.error("expected root <Session> for xml file "+xmlPath);
							return -1;
							}
						// root contains path to fasta
						if(!StringUtils.isBlank(root.getAttribute("genome"))) {
							final Path refpath = xmlPath.resolveSibling(root.getAttribute("genome"));
							if(!StringUtils.isBlank(root.getAttribute("locus"))) {
								final SAMSequenceDictionary dict=SequenceDictionaryUtils.extractRequired(refpath);
								new IntervalParser(dict).apply(root.getAttribute("locus")).ifPresent(R->roi_intervals.add(new Interval(R)));
								}
							
							
							input_files2.add(refpath);
							}
						// resources/resource/@path contains the path to the hts files
						for(Node c1=root.getFirstChild(); c1!=null;c1=c1.getNextSibling()) {
							if(c1.getNodeType()!=Node.ELEMENT_NODE) continue;
							if(!c1.getNodeName().equals("Resources")) continue;
							for(Node c2=c1.getFirstChild(); c2!=null;c2=c2.getNextSibling()) {
								if(c2.getNodeType()!=Node.ELEMENT_NODE) continue;
								final Element e2 = Element.class.cast(c2);
								if(!e2.getTagName().equals("Resource") || StringUtils.isBlank(e2.getAttribute("path"))) continue;
								input_files2.add(xmlPath.resolveSibling(e2.getAttribute("path")));
								}
							}
						}
					input_files.clear();
					input_files.addAll(input_files2);
					}
				
				for(Path htsFile:input_files) {
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
					else if(fname.endsWith(FileExtensions.COMPRESSED_VCF) || fname.endsWith(FileExtensions.COMPRESSED_VCF_BGZ)) {
						final Path index = htsFile.resolveSibling(htsFile.getFileName().toString()+ FileExtensions.TABIX_INDEX);
						if(!Files.exists(index)) {
							LOG.error("cannot find index "+index+" for vcf "+htsFile);
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
						LOG.warning("Skipping undefined file type "+htsFile);
						}
					}
				
				for(LocalFile lf1: hts_files) {
					if(!lf1.isBamOrCramOrVcf()) continue;
					final SAMSequenceDictionary dict1 = SequenceDictionaryUtils.extractRequired(lf1.path);
					lf1.reference=dict2ref.entrySet().stream()
							.filter(D->SequenceUtil.areSequenceDictionariesEqual(dict1, D.getKey()))
							.map(D->D.getValue())
							.findAny().orElse(null);
					if(lf1.reference==null) {
						LOG.error("cannot link ref for hts file "+lf1.path+" "+dict2ref+" "+dict1);
						return -1;
						}
					}
				
				
				if(this.roiBedPath!=null) {
					SAMSequenceDictionary oneDict=null;
					if(dict2ref.size()==1) {
						oneDict = dict2ref.keySet().stream().findFirst().get();
						};
					try(BedLineReader blr= BedLineReader.openPath(this.roiBedPath)) {
						blr.setContigNameConverter(oneDict);
						while(blr.hasNext()) {
							final BedLine rec =blr.next();
							final String name= rec.getOrDefault(3, "");
							roi_intervals.add(new Interval(rec.getContig(),rec.getStart(),rec.getEnd(),false,name));
							}
						}
					}
				
				// sort ROI
				if(dict2ref.size()==1) {
					Collections.sort(roi_intervals,new ContigDictComparator(dict2ref.keySet().stream().findFirst().get()).createLocatableComparator());
					};
				
				
				if(hts_files.isEmpty()) {
					LOG.error("no file was found");
					return -1;
					}
				
				server = new Server(this.port);
			
		        
		        final ServletContextHandler ctx = new ServletContextHandler();
		        ctx.setContextPath("/");
		        ctx.setResourceBase("/");
				ctx.addServlet(new ServletHolder("my", new IGVServlet(
						this.igv_js_url,
						hts_files,
						roi_intervals
						)),"/igv/*");
				
				
		        final ServletHolder files = new ServletHolder("static_file_server", DefaultServlet.class);
		        files.setInitParameter("baseResource","/");
		        files.setInitParameter("resourceBase","/");
		        files.setInitParameter("dirAllowed", "true");   // directory listing; set to false to disable
		        files.setInitParameter("welcomeServlets", "false");
		        files.setInitParameter("pathInfoOnly", "true");
		        ctx.addServlet(files, "/"+STATIC_PATH+"/*");

		        
		        
		        
		        
		        // 3. Mount it on the server
		        server.setHandler(ctx);
				
		        
				
				
				
		        LOG.info("Starting IGVServer on http://localhost:"+this.port+"/igv  . Open your browser at this URL. Press CTRL+C to stop the server.");
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
