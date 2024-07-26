/*
The MIT License (MIT)

Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.multiqc;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.google.gson.stream.JsonWriter;

/** class to generate custom content json files for multiqc */
public class MultiqcCustomContent {
	public final static String KEY_TITLE="title";
	public final static String KEY_XLAB="xlab";
	public final static String KEY_YLAB="ylab";
	public final static String KEY_SECTION_ID="section_id";
	public final static String KEY_SECTION_NAME="section_name";
	public final static String KEY_SECTION_DESC="section_description";
	public final static String KEY_PARENT_ID="parent_id";
	public final static String KEY_PARENT_NAME="parent_name";
	public final static String KEY_PARENT_DESC="parent_name_description";
	private final Map<String,String> attributes=new HashMap<>();
	
	public MultiqcCustomContent() {
		
		}
	
	/** set the attribute used for output. Key are defined as KEY** in the class */
	public MultiqcCustomContent attribute(final String key,String value) {
		this.attributes.put(key, value);
		return this;
		}
	
	protected Path validatePath(Path p) {
		String f=p.getFileName().toString();
		if(!f.endsWith("_mqc.json")) throw new IllegalArgumentException("output filename should end with _mqc.json but got "+f);
		return p;
		}
	
	protected void writeXYLab(JsonWriter w) throws IOException {
		w.name("xlab");
		w.value(this.attributes.getOrDefault(KEY_XLAB, "xlab"));
		w.name("ylab");
		w.value(this.attributes.getOrDefault(KEY_YLAB, "ylab"));
	}
	
	protected String writeSectionBegin(JsonWriter w) throws IOException {
		String section_id= this.attributes.getOrDefault(KEY_SECTION_ID,"n"+System.currentTimeMillis());
		w.name("id");
		w.value(section_id);
		w.name("section_name");
		String section_name = this.attributes.getOrDefault(KEY_SECTION_NAME, section_id);
		w.value(section_name);
		w.name("description");
		w.value(this.attributes.getOrDefault(KEY_SECTION_DESC, section_name));		
		return section_id;
		}
	
	protected void writeParentSection(JsonWriter w) throws IOException {
		final String parent_id = this.attributes.get(KEY_PARENT_ID);
		if(!StringUtils.isBlank(parent_id)) {
			w.name("parent_id");
			w.value(parent_id);
			final String parent_name =  StringUtils.ifBlank(attributes.get(KEY_PARENT_NAME),parent_id);
			w.name("parent_name");
			w.value(parent_name);
			final String parent_desc =  StringUtils.ifBlank(attributes.get(KEY_PARENT_DESC),parent_name);
			w.name("parent_description");
			w.value(parent_desc);
			}
		}
	public void writeBoxplot(Path jsonFile,Map<String,List<Number>> content) throws IOException {
		try(JsonWriter w=new JsonWriter(Files.newBufferedWriter(validatePath(jsonFile)))) {
			writeBoxplot(w,content);
			}
		}

	public void writeBoxplot(JsonWriter w,Map<String,List<Number>> content) throws IOException {
		w.beginObject();
		String section_id = writeSectionBegin(w);
		writeParentSection(w);
		
		w.name("plot_type");
		w.value("box");
		w.name("pconfig");
		w.beginObject();
		w.name("id");
		w.value("plot__id"+section_id);
		w.name("title");
		w.value(this.attributes.getOrDefault(KEY_TITLE, "title"));
		writeXYLab(w);
		
		w.endObject();
		
		w.name("data");
		w.beginObject();
		
		for(final String popName:content.keySet()) {
			w.name(popName);
			w.beginArray();
			for(Number k:content.get(popName)) {
				w.value(k);
				}
			w.endArray();
			}
			
		w.endObject();
		
		w.endObject();
		w.flush();
		}
	
	public void writeBeeswarm(Path jsonFile,Map<String,Map<String,Number>> content) throws IOException {
		try(JsonWriter w=new JsonWriter(Files.newBufferedWriter(validatePath(jsonFile)))) {
			writeBeeswarm(w,content);
			}
		}

	public void writeBeeswarm(JsonWriter w,Map<String,Map<String,Number>> content) throws IOException {
		w.beginObject();
		String section_id = writeSectionBegin(w);
		writeParentSection(w);
		
		w.name("plot_type");
		w.value("beeswarm");
		w.name("pconfig");
		w.beginObject();
		w.name("id");
		w.value("plot__id"+section_id);
		w.name("title");
		w.value(this.attributes.getOrDefault(KEY_TITLE, "title"));
		writeXYLab(w);
		
		
		w.endObject();
		
		w.name("data");
		w.beginObject();
		
		
		for(final String sn:content.keySet()) {
			w.name(sn);
			w.beginObject();
			final Map<String,Number> hash = content.get(sn);
			for(final String key2:hash.keySet()) {
				w.name(key2);
				w.value(hash.get(key2));
				}
			w.endObject();
			}
		
		w.endObject();
		
		w.endObject();
		w.flush();
		}
}
