/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.util.Collection;
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

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.samtools.SamJsonWriterFactory;
import com.google.gson.stream.JsonWriter;

public class GBrowserHtml extends AbstractGBrowserHtml
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(GBrowserHtml.class);
	
	
	public GBrowserHtml()
		{
		}

	
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		final int DEFAULT_EXTEND_INTERVAL=100;
		SamReader samReader = null;
		ZipOutputStream zout=null;
		BufferedReader bufReader = null;
		PrintWriter paramsWriter=null;
		JsonWriter paramsJsonWriter=null;
		String line;
		long snapshot_id=0L;
		try {
			final SamJsonWriterFactory samJsonWriterFactory=SamJsonWriterFactory.newInstance().
					printHeader(false).
					printAttributes(false).
					printMate(false)
					;
			if(super.getOutputFile()!=null)
				{
				if(!super.getOutputFile().getName().endsWith(".zip")) {
					return wrapException("Output file should end with *.zip");
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
			

			
			zout= new ZipOutputStream(openFileOrStdoutAsStream());
			File bamFile = null;
			File faidxFile = null;
			String sampleName = null;
			int extend_interval = 100;
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
					faidxFile= (value.isEmpty()?null:new File(value));
					}
				else if(key.equals("extend"))
					{
					extend_interval= (value.isEmpty()?DEFAULT_EXTEND_INTERVAL:Integer.parseInt(value));
					}
				else if(key.equals("position") || key.equals("interval")|| key.equals("goto"))
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
							return wrapException("bad interval :"+line);
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
					return wrapException("bad interval :"+line);
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
					samReader  =super.createSamReaderFactory().open(bamFile);
					FileWriter jsonfw = new FileWriter(tmpFile1);
					JsonWriter jsw = new JsonWriter(jsonfw);
					jsw.beginObject();
					jsw.name("reads");
					
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
						if(rec.getReadFailsVendorQualityCheckFlag()) continue;
						if(rec.getNotPrimaryAlignmentFlag()) continue;
						if(rec.getDuplicateReadFlag()) continue;
						samFileWriter.addAlignment(rec);
					}
					samRecIter.close();
					//samFileWriter.close(); NON !!
					
					jsw.endObject();
					jsw.flush();
					jsw.close();
					jsonfw.close();
					
					ZipEntry entry = new ZipEntry(super.prefix+"snapshot."+snapshot_id+".json");
					zout.putNextEntry(entry);
					IOUtils.copyTo(tmpFile1, zout);
					zout.closeEntry();
					tmpFile1.delete();
					
					paramsJsonWriter.beginObject();
					paramsJsonWriter.name("title");
					paramsJsonWriter.value(title==null?interval.toString():title);
					paramsJsonWriter.name("interval");
					paramsJsonWriter.beginObject();
						paramsJsonWriter.name("contig");
						paramsJsonWriter.value(interval.getContig());
						paramsJsonWriter.name("start");
						paramsJsonWriter.value(interval.getStart());
						paramsJsonWriter.name("end");
						paramsJsonWriter.value(interval.getEnd());
					paramsJsonWriter.endObject();
					
					if(faidxFile!=null)  {
						IndexedFastaSequenceFile faidx= new IndexedFastaSequenceFile(faidxFile);
						ReferenceSequence dna =faidx.getSubsequenceAt(interval.getContig(), interval.getStart(), interval.getEnd());
						paramsJsonWriter.name("reference");
						paramsJsonWriter.value(dna.getBaseString());
						faidx.close();
					}
					
					paramsJsonWriter.value(interval.toString());
					if(sampleName!=null && !sampleName.isEmpty()) {
						paramsJsonWriter.name("sample");
						paramsJsonWriter.value(sampleName);
					}
					paramsJsonWriter.name("href");
					paramsJsonWriter.value("snapshot."+snapshot_id+".json");
					
					paramsJsonWriter.endObject();
					}
				else if(line.toLowerCase().equals("exit") || line.toLowerCase().equals("quit"))
					{
					break;
					}
				}
			bufReader.close();bufReader=null;
			
			
			
			

			
			for(final String jssrc:new String[]{"gbrowse","hershey","samtools"}) {
				InputStream is = this.getClass().getResourceAsStream("/META-INF/js/"+jssrc+".js");
				if(is==null) {
					return wrapException("Cannot read resource /META-INF/js/"+jssrc+".js");
				}
				ZipEntry entry = new ZipEntry(super.prefix+jssrc+".js");
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
			zout.putNextEntry(new ZipEntry(super.prefix+"config.js"));
			IOUtils.copyTo(tmpFile2, zout);
			zout.closeEntry();
			tmpFile2.delete();
			
			//save index.html
			zout.putNextEntry(new ZipEntry(super.prefix+"index.html"));
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
			w.writeCharacters(getName()+":"+getVersion());
			w.writeEndElement();//title
			
			w.writeStartElement("style");
				w.writeAttribute("type", "text/css");
			w.writeEndElement();//style
			
			for(final String src:new String[]{"samtools","gbrowse","hershey","config"}) {
				w.writeStartElement("script");
				w.writeAttribute("type", "text/javascript");
				w.writeAttribute("language", "text/javascript");
				w.writeAttribute("src",src+".js");
				w.writeCharacters("");
				w.writeEndElement();//script
				}
			
			w.writeStartElement("script");
			w.writeAttribute("type", "text/javascript");
			w.writeAttribute("language", "text/javascript");
			w.writeCharacters("function init() {}");
			w.writeCharacters("window.addEventListener(\"load\",init,false);");
			w.writeEndElement();//script
			
			w.writeEndElement();//head
			
			
			w.writeStartElement("body");
			
			
			w.writeEmptyElement("canvas");
				w.writeAttribute("id", "canvasdoc");
				w.writeAttribute("width", "100");
				w.writeAttribute("height", "100");
			
			
			w.writeEndElement();//body
			w.writeEndElement();//html
			w.flush();w.close();w=null;
			
			
			zout.closeEntry();
			zout.finish();
			zout.close();
			return RETURN_OK;
		} catch(Throwable  err) {
			return wrapException(err);
		} finally
			{
			CloserUtil.close(paramsJsonWriter);
			CloserUtil.close(paramsWriter);
			CloserUtil.close(zout);
			CloserUtil.close(bufReader);
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
