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

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.samtools.SamJsonWriterFactory;
import com.google.gson.stream.JsonWriter;

/**

BEGIN_DOC




### Example



```
(...)
```






END_DOC
*/


@Program(name="forkvcf",description="Fork a VCF.")
public class GBrowserHtml extends Launcher
	{
	private static final Logger LOG = Logger.build(GBrowserHtml.class).make();

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;
	@Parameter(names={"-prefix","--prefix"},description="Zip Prefix")
	private String prefix = "gbrowse/";

	public GBrowserHtml()
		{
		}

	
	@Override
	public int doWork(List<String> args) {
		final int DEFAULT_EXTEND_INTERVAL=0;
		SamReader samReader = null;
		ZipOutputStream zout=null;
		BufferedReader bufReader = null;
		PrintWriter paramsWriter=null;
		JsonWriter paramsJsonWriter=null;
		String line;
		IndexedFastaSequenceFile faidx = null;
		long snapshot_id=0L;
		final String inputName = oneFileOrNull(args);
		try {
			final SamJsonWriterFactory samJsonWriterFactory=SamJsonWriterFactory.newInstance().
					printHeader(true).
					printAttributes(false).
					printMate(true).
					printReadQualities(false).
					closeStreamAtEnd(false)
					;
			if(this.outputFile!=null)
				{
				if(!this.outputFile.getName().endsWith(".zip")) {
					LOG.error("Output file should end with *.zip");
					return -1;
					}
				}
			
			bufReader = (inputName==null?
					IOUtils.openStreamForBufferedReader(stdin()):
					IOUtils.openURIForBufferedReading(inputName)
					);
			
			final File tmpFile1 = File.createTempFile("gbrowse.", ".tmp");
			tmpFile1.deleteOnExit();
			final File tmpFile2 = File.createTempFile("gbrowse.", ".tmp");
			tmpFile2.deleteOnExit();
			paramsWriter = new PrintWriter(tmpFile2);
			paramsWriter.print("var config=");
			paramsWriter.flush();
			paramsJsonWriter=new JsonWriter(paramsWriter);
			paramsJsonWriter.beginArray();
			

			
			zout= new ZipOutputStream(super.openFileOrStdoutAsStream(this.outputFile));
			File bamFile = null;
			File faidxFile = null;
			String sampleName = null;
			int extend_interval = DEFAULT_EXTEND_INTERVAL;
			String title=null;
			Interval interval = null;
			while( (line=bufReader.readLine())!=null) {
				LOG.info(line);
				if(line.isEmpty() || line.startsWith("#")) continue;
				final int eq=line.indexOf("=");
				final String key =   (eq==-1?"":line.substring(0, eq).toLowerCase().trim());
				final String value = (eq==-1?"":line.substring(eq+1).trim());
				
				if(key.equals("bam"))
					{
					if(samReader!=null) samReader.close();
					samReader=null;
					bamFile= (value.isEmpty()?null:new File(value));
					}
				else if(key.equals("sample"))
					{
					sampleName= (value.isEmpty()?null:value);
					}
				else if(key.equals("title"))
					{
					title= (value.isEmpty()?null:value);
					}
				else if(key.equals("ref") || key.equals("fasta"))
					{
					if( faidx != null) faidx.close();
					faidx=null;
					faidxFile= (value.isEmpty()?null:new File(value));
					}
				else if(key.equals("extend"))
					{
					extend_interval= (value.isEmpty()?DEFAULT_EXTEND_INTERVAL:Integer.parseInt(value));
					}
				else if(key.equals("position") || key.equals("location") || key.equals("interval")|| key.equals("goto"))
					{
					Pattern pat1 = Pattern.compile("([^\\:]+)\\:([\\d,]+)");
					Matcher matcher = pat1.matcher(value);
					if(matcher.matches())
						{
						String c= matcher.group(1);
						int pos = Integer.parseInt(matcher.group(2).replaceAll("[,]",""));
						interval = new Interval(c,
								Math.max(1, pos-extend_interval),
								pos+extend_interval
								);
						continue;
						}
					pat1 = Pattern.compile("([^\\:]+)\\:([\\d,]+)\\-([\\d,]+)");
					matcher = pat1.matcher(value);
					if(matcher.matches())
						{
						String c= matcher.group(1);
						int B = Integer.parseInt(matcher.group(2).replaceAll("[,]",""));
						int E = Integer.parseInt(matcher.group(3).replaceAll("[,]",""));
						if(B>E) {
							LOG.error("bad interval :"+line);
							return -1;
						}
						interval = new Interval(c,Math.max(1, B-extend_interval),E+extend_interval);
						continue;
						}
					pat1 = Pattern.compile("([^\\:]+)\\:([\\d,]+)\\+([\\d,]+)");
					matcher = pat1.matcher(value);
					if(matcher.matches())
						{
						String c= matcher.group(1);
						int B = Integer.parseInt(matcher.group(2).replaceAll("[,]",""));
						int x = Integer.parseInt(matcher.group(3).replaceAll("[,]",""));
						interval = new Interval(c,Math.max(1, B-(x+extend_interval)),B+(x+extend_interval));
						continue;
						}
					LOG.error("bad interval :"+line);
					return -1;
					}
				else if(line.toLowerCase().equals("snapshot"))
					{
					if(interval==null) {
						LOG.error("No interval defined!");
						continue;
						}
					if(bamFile==null) {
						LOG.error("No BAM file defined!");
						continue;
						}
					++snapshot_id;
					
					LOG.info("open samFile "+bamFile);
					if(samReader==null) 
						{
						samReader  =super.createSamReaderFactory().open(bamFile);
						}
					FileWriter jsonFileWriter = new FileWriter(tmpFile1);
					JsonWriter jsw = new JsonWriter(jsonFileWriter);
					jsw.beginObject();
					jsw.name("interval");
					jsw.beginObject();
						jsw.name("contig");
						jsw.value(interval.getContig());
						jsw.name("start");
						jsw.value(interval.getStart());
						jsw.name("end");
						jsw.value(interval.getEnd());
					jsw.endObject();

					if(faidxFile!=null)  {
						if( faidx == null) {
							faidx= new IndexedFastaSequenceFile(faidxFile);
							}
						ReferenceSequence dna =faidx.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
						jsw.name("reference");
						jsw.value(dna.getBaseString());
						
					}


					jsw.name("sam");
					
					SAMFileWriter samFileWriter=samJsonWriterFactory.open(samReader.getFileHeader(), jsw);				
					SAMRecordIterator samRecIter = samReader.queryOverlapping(interval.getContig(), interval.getStart(), interval.getEnd());
					while(samRecIter.hasNext()) {
						final SAMRecord rec= samRecIter.next();
						if(sampleName!=null && !sampleName.isEmpty())
							{
							SAMReadGroupRecord srg = rec.getReadGroup();
							if(srg==null) continue;
							if(!sampleName.equals(srg.getSample())) continue;
							}
						if(rec.getReadUnmappedFlag()) continue;
						
						samFileWriter.addAlignment(rec);
					}
					samRecIter.close();
					samFileWriter.close();
					
					jsw.endObject();
					jsw.flush();
					jsw.close();
					jsonFileWriter.close();
					
					ZipEntry entry = new ZipEntry(this.prefix+"_snapshot."+String.format("%05d",snapshot_id)+".json");
					zout.putNextEntry(entry);
					IOUtils.copyTo(tmpFile1, zout);
					zout.closeEntry();
					tmpFile1.delete();
					
					paramsJsonWriter.beginObject();
					paramsJsonWriter.name("title");
					paramsJsonWriter.value(title==null?
						interval.getContig()+":"+interval.getStart()+"-"+interval.getEnd()+(sampleName==null?"":" "+sampleName)
						:title);
					paramsJsonWriter.name("interval");
					paramsJsonWriter.beginObject();
						paramsJsonWriter.name("contig");
						paramsJsonWriter.value(interval.getContig());
						paramsJsonWriter.name("start");
						paramsJsonWriter.value(interval.getStart());
						paramsJsonWriter.name("end");
						paramsJsonWriter.value(interval.getEnd());
					paramsJsonWriter.endObject();
					

					if(sampleName!=null && !sampleName.isEmpty()) {
						paramsJsonWriter.name("sample");
						paramsJsonWriter.value(sampleName);
					}
					paramsJsonWriter.name("href");
					paramsJsonWriter.value("_snapshot."+String.format("%05d",snapshot_id)+".json");
					
					paramsJsonWriter.endObject();
					}
				else if(line.toLowerCase().equals("exit") || line.toLowerCase().equals("quit"))
					{
					break;
					}
				}
			bufReader.close();bufReader=null;
			
			
			
			

			
			for(final String jssrc:new String[]{"gbrowse","hershey","samtools","com.github.lindenb.jvarkit.tools.misc.GBrowserHtml"}) {
				InputStream is = this.getClass().getResourceAsStream("/META-INF/js/"+jssrc+".js");
				if(is==null) {
					LOG.error("Cannot read resource /META-INF/js/"+jssrc+".js");
					return -1;
				}
				ZipEntry entry = new ZipEntry(this.prefix+jssrc+".js");
				entry.setComment("JAVASCRIPT SOURCE "+jssrc);
				zout.putNextEntry(entry);
				IOUtils.copyTo(is, zout);
				CloserUtil.close(is);
				zout.closeEntry();				
			}
			
			
			//save params.js
			paramsJsonWriter.endArray();
			paramsJsonWriter.flush();
			paramsWriter.println(";");
			paramsWriter.flush();
			paramsJsonWriter.close();
			paramsWriter.close();
			zout.putNextEntry(new ZipEntry(this.prefix+"config.js"));
			IOUtils.copyTo(tmpFile2, zout);
			zout.closeEntry();
			tmpFile2.delete();
			
			//save index.html
			zout.putNextEntry(new ZipEntry(this.prefix+"index.html"));
			XMLOutputFactory xof = XMLOutputFactory.newFactory();
			XMLStreamWriter w  = xof.createXMLStreamWriter(zout,"UTF-8");
			w.writeStartElement("html");
			w.writeStartElement("head");
			w.writeEmptyElement("meta");
				w.writeAttribute("http-equiv", "Content-Type");
				w.writeAttribute("content", "text/html; charset=utf-8");
			w.writeEmptyElement("meta");
				w.writeAttribute("http-equiv", "author");
				w.writeAttribute("content", "Pierre Lindenbaum Phd @yokofakun");
			
			w.writeStartElement("title");
			w.writeCharacters(getProgramName()+":"+getVersion());
			w.writeEndElement();//title
			
			w.writeStartElement("style");
				w.writeAttribute("type", "text/css");
			
			w.writeCharacters("body	{ color:rgb(50,50,50); margin:20px; padding:20px; font: 12pt Arial, Helvetica, sans-serif; }\n");
			w.writeCharacters("label { text-align:right;	 }\n");
			w.writeCharacters("button	{ border: 1px solid; background-image:-moz-linear-gradient( top, gray, lightgray ); }\n");
			w.writeCharacters("canvas { image-rendering:auto;}\n");
			w.writeCharacters(".me 	{ padding-top:100px; font-size:80%; }\n");
			w.writeEndElement();//style
			
			for(final String src:new String[]{"samtools","gbrowse","hershey","config","com.github.lindenb.jvarkit.tools.misc.GBrowserHtml"}) {
				w.writeStartElement("script");
				w.writeAttribute("type", "text/javascript");
				w.writeAttribute("language", "text/javascript");
				w.writeAttribute("src",src+".js");
				w.writeCharacters("");
				w.writeEndElement();//script
				}
			
			
			w.writeEndElement();//head
			
			
			w.writeStartElement("body");
			
			w.writeStartElement("div");
			w.writeStartElement("button");
				w.writeAttribute("onclick", "changemenu(-1)");
				w.writeCharacters("prev");
			w.writeEndElement();//button
			
			w.writeStartElement("select");
			w.writeAttribute("id", "menu");
			w.writeEndElement();//menu
			
			w.writeStartElement("button");
				w.writeAttribute("onclick", "changemenu(+1)");
				w.writeCharacters("next");
			w.writeEndElement();//button
			w.writeEndElement();//div
			
			w.writeStartElement("div");
			w.writeAttribute("id", "flags");
			w.writeEndElement();//div
			
			
			w.writeStartElement("div");
			w.writeAttribute("style", "text-align:center;");
			
			w.writeStartElement("div");
			w.writeAttribute("style", "font-size:200%;margin:10px;");
			w.writeAttribute("id", "browserTitle");
			w.writeEndElement();//div
			
			w.writeStartElement("div");
			w.writeEmptyElement("canvas");
				w.writeAttribute("id", "canvasdoc");
				w.writeAttribute("width", "100");
				w.writeAttribute("height", "100");
			w.writeEndElement();//div
			w.writeEndElement();//div
			
			w.writeEmptyElement("hr");
			w.writeStartElement("div");
			w.writeAttribute("class", "me");
			w.writeCharacters("Pierre Lindenbaum PhD. ");
			w.writeStartElement("a");
			w.writeAttribute("href", "https://github.com/lindenb/jvarkit");
			w.writeCharacters("https://github.com/lindenb/jvarkit");
			w.writeEndElement();
			w.writeCharacters(". Tested with Firefox 45.0");
			w.writeEndElement();
			
			
			w.writeEndElement();//body
			w.writeEndElement();//html
			w.flush();w.close();w=null;
			
			
			zout.closeEntry();
			zout.finish();
			zout.close();
			return RETURN_OK;
		} catch(Throwable  err) {
			LOG.error(err);
			return -1;
		} finally
			{
			CloserUtil.close(paramsJsonWriter);
			CloserUtil.close(paramsWriter);
			CloserUtil.close(zout);
			CloserUtil.close(bufReader);
			CloserUtil.close(samReader);
			CloserUtil.close(faidx);
			}
		}
	
	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new GBrowserHtml().instanceMainWithExit(args);
		}

}
