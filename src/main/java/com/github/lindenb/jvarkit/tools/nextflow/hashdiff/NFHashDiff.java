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
package com.github.lindenb.jvarkit.tools.nextflow.hashdiff;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonWriter;
/**
BEGIN_DOC

Nextflow run **must** be launcher with `-dump-hashes json`

END_DOC
*/
@Program(name="nfhashdiff",
description="Compare .nextflow.log",
keywords={"json"},
creationDate="20251119",
modificationDate="20251119",
jvarkit_amalgamion = false,
jvarkit_hidden = true,
menu="Utilities"
)
public class NFHashDiff extends Launcher {
	private static final Logger LOG = Logger.of(NFHashDiff.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile = null;
	@Parameter(names={"-1"},description="hide unique to file1")
	private boolean hide_file1=false;
	@Parameter(names={"-2"},description="hide unique to file2")
	private boolean hide_file2=false;
	@Parameter(names={"-t"},description="look only at hashes value that value")
	private String search_value=null;
	@Parameter(names={"-g","--grep"},description="Search for all those strings")
	private List<String> grep_items = new ArrayList<>();


	private int compare(JsonArray a,JsonArray b) {
		int i = Integer.compare(a.size(), b.size());
		if(i!=0) return i;
		for(int x=0;x< a.size();++x) {
			JsonObject oa = a.get(x).getAsJsonObject();
			String sa = oa.get("hash").getAsString();
			JsonObject ob = b.get(x).getAsJsonObject();
			String sb = ob.get("hash").getAsString();
			i = sa.compareTo(sb);
			if(i!=0) return i;
			}
		return 0;
		}
	private boolean equals(JsonArray a,JsonArray b) {
		return compare(a,b)==0;
		}
	
	private List<JsonArray> readLog(Path p) throws IOException {
		try(BufferedReader br = IOUtils.openPathForBufferedReading(p) ) {
			final List<JsonArray> elements = new ArrayList<>();
			for(;;) {
				String line=br.readLine();
				if(line==null) break;
				if(!line.contains("cache hash:")) continue;
				if(!line.endsWith("[")) {
					LOG.warn("line doesn't end with '[' . Did you use '-dump-hashes json' ? "+line); 
					continue;
					}
				final StringBuilder sb=new StringBuilder("[");
				for(;;) {
					line=br.readLine();
					if(line==null) {
						throw new IOException("unclosed json");
						}
					sb.append(line);
					if(line.equals("]")) {
						JsonParser r=new JsonParser();
						JsonArray array= r.parse(sb.toString()).getAsJsonArray();
						boolean ok=true;
						
						if(!StringUtils.isBlank(this.search_value)) {
							ok=false;
							for(int x=0;x< array.size();++x) {
								JsonObject oa = array.get(x).getAsJsonObject();
								String sa = oa.get("value").getAsString();
								if(sa.equals(this.search_value)) {
									ok=true;
									break;
									}
								}
							}
						if(ok && !grep_items.isEmpty()) {
							int n=0;
							for(String s:this.grep_items) {
								boolean found=false;
								for(int x=0;x< array.size();++x) {
									final JsonObject oa = array.get(x).getAsJsonObject();
									String sa = oa.get("value").getAsString();
									if(sa.contains(s)) {
										found=true;
										break;
										}
									sa = oa.get("hash").getAsString();
									if(sa.contains(s)) {
										found=true;
										break;
										}
									}
								n+=(found?1:0);
								if(!found) break;
								}
							ok = (n==this.grep_items.size());
							}
						
						if(ok) {
							elements.add(array);
							}
						break;
						}
					}
				}
			LOG.info("in "+p+" got "+elements.size()+" element(s).");
			if(elements.isEmpty()) {
				LOG.warn("no element found. Did you use '-dump-hashes json'  ?");
				}
			Collections.sort(elements,(A,B)->compare(A,B));
			return elements;
			}
		}
	private void print(JsonWriter w,Gson gson, List<JsonArray> L) throws IOException {
		w.beginArray();
		for(int i=0;i< L.size();i++) {
			w.jsonValue(gson.toJson( L.get(i)));
			}
		w.endArray();
		}
	@Override
	public int doWork(List<String> args) {
		try {
			if(args.size()!=1 && args.size()!=2) {
				LOG.error("exepected one or two files");
				return -1;
				}
			final List<JsonArray> L1 = readLog(Paths.get(args.get(0)));
			
			
			if(args.size()==1) {
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outFile)) {
					Gson gson = new GsonBuilder().setPrettyPrinting().create();
					try(JsonWriter w=gson.newJsonWriter(out)) {
						print(w,gson,L1);
						}
					}
				}
			else
				{
				final List<JsonArray> L2 = readLog(Paths.get(args.get(1)));
				
				while(!L1.isEmpty() && !L2.isEmpty()) {
					final JsonArray oa = L1.get(0);
					final JsonArray ob = L2.get(0);
					if(equals(oa,ob)) {
						L1.remove(0);
						L2.remove(0);
						}
					else
						{
						break;
						}
					}
				
				int x=0;
				while(x< L1.size()) {
					JsonArray oa = L1.get(x);
					int y=0;
					while(y < L2.size()) {
						JsonArray ob = L2.get(y);
						if(equals(oa,ob)) break;
						y++;
						}
					if(y< L2.size()) {
						L1.remove(x);
						L2.remove(y);
						}
					else
						{
						x++;
						}
					}
				
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outFile)) {
					Gson gson = new GsonBuilder().setPrettyPrinting().create();
					try(JsonWriter w=gson.newJsonWriter(out)) {
						w.beginArray();
						
						
						
						
						if(!hide_file1) {
							w.beginObject();
							w.name("description");
							w.value("Unique to "+args.get(0));
							w.name("file");
							w.value(args.get(0));
							w.name("items");
							print(w,gson,L1);
							w.endObject();
							}
						if(!hide_file2) {
							w.beginObject();
							w.name("description");
							w.value("Unique to "+args.get(1));
							w.name("file");
							w.value(args.get(1));
							w.name("items");
							print(w,gson,L2);
							w.endObject();
							}
						w.endArray();
						w.flush();
					}
				}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new NFHashDiff().instanceMainWithExit(args);
	}


}
