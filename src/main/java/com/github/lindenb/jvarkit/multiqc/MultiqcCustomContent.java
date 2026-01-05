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
package com.github.lindenb.jvarkit.multiqc;

import java.awt.geom.Point2D;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.Maps;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
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
	
	private String getChartTitle() {
		String s= this.attributes.getOrDefault(KEY_TITLE,"");
		if(StringUtils.isBlank(s)) s=this.attributes.getOrDefault(KEY_SECTION_NAME,"");
		if(StringUtils.isBlank(s)) s=this.attributes.getOrDefault(KEY_SECTION_ID,"");
		if(StringUtils.isBlank(s)) s="title";
		return s;
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
	
	protected JsonObject putXYLab(JsonObject o) throws IOException {
		o.addProperty("xlab",this.attributes.getOrDefault(KEY_XLAB, "xlab"));
		o.addProperty("ylab",this.attributes.getOrDefault(KEY_YLAB, "ylab"));
		return o;
		}
	
	protected void writeXYLab(JsonWriter w) throws IOException {
		w.name("xlab");
		w.value(this.attributes.getOrDefault(KEY_XLAB, "xlab"));
		w.name("ylab");
		w.value(this.attributes.getOrDefault(KEY_YLAB, "ylab"));
	}
	
	
	protected JsonObject putSectionBegin(JsonObject o) throws IOException {
		String section_id= this.attributes.getOrDefault(KEY_SECTION_ID,"n"+System.currentTimeMillis());
		o.addProperty("id",section_id);
		
		String section_name = this.attributes.getOrDefault(KEY_SECTION_NAME, section_id);
		o.addProperty("section_name",section_name);
		o.addProperty("description",this.attributes.getOrDefault(KEY_SECTION_DESC, section_name));		
		return o;
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
			w.setIndent("    ");
			w.setLenient(true);
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
		w.value(getChartTitle());
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
	public void writeBarGraph(JsonWriter w,Map<String,Number> hash) throws IOException {
		writeBarGraph(w,"Count",false,hash);
		}
	/**
	packed: put everything in the same bar, otherwise, one bar for each key
	*/
	
	public void writeBarGraph(JsonWriter w,final String countTitle,boolean packed, Map<String,Number> hash) throws IOException {
		final Map<String,Map<String,Number>> hash2=new LinkedHashMap<>();
		if(!packed) {
			for(String key:hash.keySet()) {
				hash2.put(key, Maps.of(countTitle, hash.get(key)));
				}
			}
		else
			{
			hash2.put(countTitle, hash);
			}
		writeMultiBarGraph(w,hash2);
		}
	
	public void writeMultiBarGraph(JsonWriter w,Map<String,Map<String,Number>> hash1) throws IOException {
		w.beginObject();
		final String section_id = writeSectionBegin(w);
		writeParentSection(w);
		
		w.name("plot_type");
		w.value("bargraph");
		w.name("pconfig");
		w.beginObject();
			w.name("id");
			w.value("plot__id"+section_id);
			w.name("title");
			w.value(getChartTitle());
			writeXYLab(w);
		w.endObject();
		
		w.name("data");
		
		w.beginObject();
		for(final String key1:hash1.keySet()) {
			final Map<String,Number> hash2 = hash1.get(key1);
			w.name(key1);
			
			w.beginObject();
			for(final String key2:hash2.keySet()) {
				w.name(key2);
				w.value(hash2.get(key2));
				}
			w.endObject();
			}
		w.endObject();
		
		
		w.endObject();
		w.flush();
		}
	
	public void writeXY(JsonWriter w,List<Point2D> points) throws IOException {
		int n=0;
		final Map<String,Point2D> hash2point=new LinkedHashMap<>(points.size());
		for(Point2D pt:points) {
			hash2point.put(String.valueOf(n++), pt);
			}
		writeXY(w,hash2point);
		}
	
	public JsonObject makeXY(final Map<String,Point2D> hash2point) throws IOException {
	
		final JsonObject pconfig = new JsonObject();
		pconfig.addProperty("id", "plot_id"+"TODO");
		pconfig.addProperty("title", getChartTitle());
		
		final JsonObject o1 = new JsonObject();
		
		o1.addProperty("plot_type", "scatter");
		o1.add("pconfig", pconfig);
		
		
		final JsonObject data = new JsonObject();
		o1.add("data", data);
		
		for(final String key1:hash2point.keySet()) {
			final Point2D pt = hash2point.get(key1);
			final JsonObject pto = new JsonObject();
			pto.addProperty("x", pt.getX());
			pto.addProperty("y", pt.getY());
			data.add(key1, pto);
			}
		
		return o1;
		}
	
	public void writeXY(JsonWriter w,Map<String,Point2D> hash2point) throws IOException {
		w.beginObject();
		final String section_id = writeSectionBegin(w);
		writeParentSection(w);
		
		w.name("plot_type");
		w.value("scatter");
		w.name("pconfig");
		w.beginObject();
			w.name("id");
			w.value("plot__id"+section_id);
			w.name("title");
			w.value(getChartTitle());
			writeXYLab(w);
		w.endObject();
		
		w.name("data");
		
		w.beginObject();
		for(final String key1:hash2point.keySet()) {
			final Point2D pt = hash2point.get(key1);
			w.name(key1);
			w.beginObject();
			
				w.name("x");
				w.value(pt.getX());
				w.name("y");
				w.value(pt.getY());

			w.endObject();
			}
		w.endObject();
		
		
		w.endObject();
		w.flush();
		}
	
	public void writeLineGraph(JsonWriter w,List<Point2D> points) throws IOException {
		writeLineGraph(w,"points",points);
		}
	
	public void writeLineGraph(JsonWriter w,String name,List<Point2D> points) throws IOException {
		writeLineGraph(w, Maps.of(name,points));
		}
	
	public JsonObject toLineGraph(Map<String,List<Point2D>> hash2point) {
		String section_id="";
		JsonObject pconfig = new JsonObject();
		pconfig.addProperty("id","plot__id"+section_id);
		
		JsonObject data = new JsonObject();
		for(final String key1:hash2point.keySet()) {
			JsonObject ptobj = new JsonObject();
			for(final Point2D pt : hash2point.get(key1)) {
				ptobj.addProperty(String.valueOf(pt.getX()), pt.getY());
				}
			data.add(key1,ptobj);
			}
		
		final JsonObject o1 = new JsonObject();
		o1.addProperty("plot_type", "linegraph");
		o1.add("pconfig",pconfig);
		o1.add("data",data);
		return o1;
		}
	
	public void writeLineGraph(JsonWriter w,Map<String,List<Point2D>> hash2point) throws IOException {
		w.beginObject();
		final String section_id = writeSectionBegin(w);
		writeParentSection(w);
		
		w.name("plot_type");
		w.value("linegraph");
		w.name("pconfig");
		w.beginObject();
			w.name("id");
			w.value("plot__id"+section_id);
			w.name("title");
			w.value(getChartTitle());
			writeXYLab(w);
		w.endObject();
		
		w.name("data");
		
		w.beginObject();
		for(final String key1:hash2point.keySet()) {
			w.name(key1);
			w.beginObject();
			for(final Point2D pt : hash2point.get(key1)) {
				w.name(String.valueOf(pt.getX()));
				w.value(pt.getY());
				}
			w.endObject();
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

	public JsonElement toBeeswarm(final Map<String,Map<String,Number>> content) {
		final JsonObject data = new JsonObject();
		for(final String sn:content.keySet()) {
			final Map<String,Number> hash = content.get(sn);
			final JsonObject h2=new JsonObject();
			
			for(final String key2:hash.keySet()) {
				h2.addProperty(key2,hash.get(key2));
				}
		
			data.add(sn,h2);
			}
		
		final JsonObject chart = new JsonObject();
		chart.add("data", data);
		
		
		return chart;
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
		w.value(getChartTitle());
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
