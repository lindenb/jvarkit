/*
The MIT License (MIT)

Copyright (c) 2023 Pierre Lindenbaum

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

*/package com.github.lindenb.jvarkit.tools.samplesrdf;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.jena.rdf.model.Literal;
import org.apache.jena.rdf.model.Model;
import org.apache.jena.rdf.model.ModelFactory;
import org.apache.jena.rdf.model.Property;
import org.apache.jena.rdf.model.RDFNode;
import org.apache.jena.rdf.model.Resource;
import org.apache.jena.rdf.model.ResourceFactory;
import org.apache.jena.vocabulary.RDF;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.stream.CollectorsUtils;
import com.github.lindenb.jvarkit.stream.HtsCollectors;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;

import htsjdk.samtools.util.StringUtil;

public class SamplesRDF extends Launcher {
	
	private final Map<String,Property> id2property = new HashMap<>();
	private static final String NS = "";
	private static final Resource TYPE_SAMPLE = ResourceFactory.createResource(NS+"Sample");
	private static final Resource TYPE_GROUP = ResourceFactory.createResource(NS+"Group");
	private static final Property PROP_ID = ResourceFactory.createProperty(NS,"id");
	
	
	private final Model model = ModelFactory.createDefaultModel();
	
	
	private Property toProperty(String s) {
		final Property p = this.id2property.get(s.toLowerCase());
		if(p==null) throw new IllegalArgumentException("Undefined property \""+s+"\". Available are "+ String.join(", ", id2property.keySet()));
		return p;
		}
	
	private RDFNode toValue(Property p,String s) {
		if(p.equals(RDF.type)) {
			if(s.equalsIgnoreCase("sample")) {
				return TYPE_SAMPLE;
				}
			else if(s.equalsIgnoreCase("group")) {
				return TYPE_GROUP;
				}
			throw new IllegalArgumentException("unsupported "+p+" : "+s);
 			}
		if(p.equals(PROP_ID)) {
			s=s.toUpperCase();
			if(!s.matches("[A-Za-z0-9_]+")) throw new IllegalArgumentException("not valid pattern for "+p+"="+s);
			}
		return ResourceFactory.createPlainLiteral(s);
		}
	
	private Map.Entry<Property, RDFNode> split(final String s) {
		int i= s.indexOf(":");
		if(i==-1) i=s.indexOf("=");
		if(i==-1) throw new IllegalArgumentException("Cannot find delimiter in " + s);
		final String key= s.substring(0,i).trim().toLowerCase();
		if(key.isEmpty()) throw new IllegalArgumentException("Empty key in " + s);
		final String value = s.substring(i+1).trim();
		if(value.isEmpty()) throw new IllegalArgumentException("Empty value in " + s);
		final Property p=toProperty(key);
		return new AbstractMap.SimpleEntry<>(p, toValue(p,value));
		}
	
	private Resource toResource(RDFNode r) {
		return r.asResource();
		}
	private Literal toLiteral(RDFNode r) {
		return r.asLiteral();
		}
	
	private void scan(BufferedReader br) throws IOException {
		final List<Map.Entry<Property, RDFNode>> props = new ArrayList<>();
		for(;;) {
			String line = br.readLine();
			if(StringUtil.isBlank(line)) {
				if(!props.isEmpty()) {
					
					final Literal id = 
							props.stream().
							filter(P->P.getKey().equals(PROP_ID)).
							map(KV->toLiteral(KV.getValue())).
							collect(CollectorsUtils.one())
							;
					final Resource subject = this.model.createResource(NS+"#"+id.getString());
					if(this.model.contains(null, PROP_ID, id)) {
						throw new IOException("Duplicate definition of entity with id "+id);
						}
					if(this.model.containsResource(subject)) {
						throw new IOException("Duplicate definition of subject "+subject);
						}
					for(Map.Entry<Property, RDFNode> prop:props) {
						this.model.add(subject,prop.getKey(), prop.getValue());
						}
					}
				if(line==null) break;
				props.clear();
				}
			if(line.startsWith("#")) continue;
			props.add(split(line));
			}
		}
	
	@Override
	public int doWork(List<String> args) {
		try {
			this.id2property.put("type", RDF.type);
			this.id2property.put("id", PROP_ID);
			
			final List<Path> paths = IOUtils.unrollPaths(args);
			if(paths.isEmpty()) {
				try(BufferedReader br = IOUtils.openStreamForBufferedReader(stdin())) {
					scan(br);
					}
				}
			else
				{
				for(Path path: paths) {
					try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
						scan(br);
						}
					}
				}
			return 0;
			}
		catch(Throwable err) {
			return -1;
			}
		finally
			{
			this.model.close();
			}
		}
	
	public static void main(String[] args) {
		new SamplesRDF().instanceMainWithExit(args);

	}

}
