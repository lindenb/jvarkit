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
package com.github.lindenb.jvarkit.tools.nextflow.nfhashtools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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

# Motivation

this tool looks  for the log file of nextflow ( `.nextflow.log') and generate a json output.

Requirement: Nextflow run **must** be launched with `-dump-hashes json`

If there is no file, the input is read on stdin.

If there is only one file, the hashes are collected, filtered and printed as JSON.

If there are two files. Both files are read and filtered. The JSON output contains
- a list of objects uniques to the file1, 
- a list of objects uniques to the file2
- a list of common objects common to file1 and file2

If there only one object remaining in each list, an list of the common/discordant fields is also printed.


# Example


```
$ sdiff newflow.1.log newflow.2.log
Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo   Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo
    {                                                               {
        "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",           |         "hash": "ebc58ab2cb4848d04ec23d83f7ddf985",
        "type": "java.util.UUID",                                       "type": "java.util.UUID",
        "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"                 "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
    },                                                              },
    {                                                               {
        "hash": "f88ed0fa50638ac1efa458e70e21ebcc",                     "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
        "type": "java.lang.String",                                     "type": "java.lang.String",
        "value": "Hello"                                      |         "value": "Bonjour"
    }                                                               }
]                                                               ]
Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo   Nov-19 22:23:53.071 [Actor Thread 96] INFO  nextflow.processo
    {                                                               {
        "hash": "1c72b78fd00c568c524095ff174b1a6f",                     "hash": "1c72b78fd00c568c524095ff174b1a6f",
        "type": "java.util.UUID",                                       "type": "java.util.UUID",
        "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"                 "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
    },                                                              },
    {                                                               {
        "hash": "cb2d18e15fc931db472229e671a76351",                     "hash": "cb2d18e15fc931db472229e671a76351",
        "type": "java.lang.String",                                     "type": "java.lang.String",
        "value": "World"                                                "value": "World"
    }                                                               }
]                                                               ]

```

one file:

```
$ java -jar jvarkit.jar nfhashtools nextflow.1.log 2> /dev/null 
[
  [
  {
    "hash": "1c72b78fd00c568c524095ff174b1a6f",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "cb2d18e15fc931db472229e671a76351",
    "type": "java.lang.String",
    "value": "World"
  }
],
  [
  {
    "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Hello"
  }
]
]
```

one file with a filter:
```
$ java -jar jvarkit.jar nfhashtools -t Hello nextflow.1.log 2> /dev/null 
[
  [
  {
    "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Hello"
  }
]
]
```

with two files 

```
$ java -jar jvarkit.jar nfhashtools  nextflow.1.log nextflow.2.log 2> /dev/null 
[
  {
    "name": "L1",
    "description": "Unique to nextflow.1.log",
    "file": "nextflow.1.log",
    "items": [
      [
  {
    "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Hello"
  }
]
    ]
  },
  {
    "name": "L2",
    "description": "Unique to nextflow.2.log",
    "file": "nextflow.2.log",
    "items": [
      [
  {
    "hash": "ebc58ab2cb4848d04ec23d83f7ddf985",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
    "type": "java.lang.String",
    "value": "Bonjour"
  }
]
    ]
  },
  {
    "name": "L3",
    "description": "common to nextflow.1.log and nextflow.2.log",
    "file": [
      "nextflow.1.log,
      "nextflow.2.log"
    ],
    "items": [
      [
  {
    "hash": "1c72b78fd00c568c524095ff174b1a6f",
    "type": "java.util.UUID",
    "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
  },
  {
    "hash": "cb2d18e15fc931db472229e671a76351",
    "type": "java.lang.String",
    "value": "World"
  }
]
    ]
  },
  {
    "name": "uniq",
    "description": "common/discordant values shared between the unique item of L1 and the unique item of L2",
    "common": [],
    "uniq_1": [
      {
  "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
  "type": "java.lang.String",
  "value": "Hello"
},
      {
  "hash": "5350de1bf7fd1f286a1ad5baa7af30f5",
  "type": "java.util.UUID",
  "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
}
    ],
    "uniq_2": [
      {
  "hash": "f88ed0fa50638ac1efa458e70e21ebcc",
  "type": "java.lang.String",
  "value": "Bonjour"
},
      {
  "hash": "ebc58ab2cb4848d04ec23d83f7ddf985",
  "type": "java.util.UUID",
  "value": "7bae3014-6912-4d03-8e91-0db03ce8556b"
}
    ]
  }
]
```


END_DOC
*/
@Program(name="nfhashtools",
description="A tool to print and compare hashes .nextflow.log files",
keywords={"json","nf","nextflow"},
creationDate="20251119",
modificationDate="20251120",
jvarkit_amalgamion = true,
jvarkit_hidden = false,
generate_doc = true,
menu="Utilities"
)
public class NFHashTools extends Launcher {
	private static final Logger LOG = Logger.of(NFHashTools.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outFile = null;
	@Parameter(names={"-1"},description="hide items unique to file1")
	private boolean hide_file1=false;
	@Parameter(names={"-2"},description="hide items unique to file2")
	private boolean hide_file2=false;
	@Parameter(names={"-3"},description="hide common items in  file1 and file2")
	private boolean hide_file3=false;
	@Parameter(names={"-u"},description="hide common fields of the unique item in  file1 and file2")
	private boolean hide_unique =false;
	@Parameter(names={"-t"},description="look only at 'value' matching that string")
	private String search_value=null;
	@Parameter(names={"-g","--grep"},description="Search for all those (sub)strings in 'hash' and 'value'. Case sensible . The item matching all the regex are kept.")
	private List<String> grep_items = new ArrayList<>();


	private int compare(final JsonArray a,final JsonArray b) {
		int i = Integer.compare(a.size(), b.size());
		if(i!=0) return i;
		for(int x=0;x< a.size();++x) {
			final JsonObject oa = a.get(x).getAsJsonObject();
			final String sa = oa.get("hash").getAsString();
			final JsonObject ob = b.get(x).getAsJsonObject();
			final String sb = ob.get("hash").getAsString();
			i = sa.compareTo(sb);
			if(i!=0) return i;
			}
		return 0;
		}
	private boolean equals(final JsonArray a,final JsonArray b) {
		return compare(a,b)==0;
		}
	private List<JsonArray> readLog(final Path p) throws IOException {
		try(BufferedReader br = IOUtils.openPathForBufferedReading(p) ) {
			return readLog(br,p==null?"-":p.toString());
			}
		}
		
	private List<JsonArray> readLog(final BufferedReader br,final String filename) throws IOException {
			boolean got_one=false;
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
						final JsonParser r=new JsonParser();
						final JsonArray array= r.parse(sb.toString()).getAsJsonArray();
						got_one=true;
						boolean ok=true;
						
						if(!StringUtils.isBlank(this.search_value)) {
							ok=false;
							for(int x=0;x< array.size();++x) {
								final JsonObject oa = array.get(x).getAsJsonObject();
								final String sa = oa.get("value").getAsString();
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
			LOG.info("in "+(filename==null?"-":filename)+" got "+elements.size()+" element(s).");
			if(!got_one) {
				LOG.warn("no element found. Did you use '-dump-hashes json'  ?");
				}
			Collections.sort(elements,(A,B)->compare(A,B));
			return elements;
			}
	private void print(final JsonWriter w,final Gson gson, final List<JsonArray> L) throws IOException {
		w.beginArray();
		for(int i=0;i< L.size();i++) {
			w.jsonValue(gson.toJson( L.get(i)));
			}
		w.endArray();
		}
	@Override
	public int doWork(final List<String> args) {
		try {
			if(args.size()>2) {
				LOG.error("exepected zero, one or two .nextflow.log files");
				return -1;
				}
			final List<JsonArray> L1 = args.isEmpty() || args.get(0).equals("-")?
					readLog(IOUtils.openStreamForBufferedReader(stdin()),"<stdin>"):
					readLog(Paths.get(args.get(0)));
			
			final Gson gson = new GsonBuilder().setPrettyPrinting().create();

			if(args.size()<2) {
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outFile)) {
					try(JsonWriter w=gson.newJsonWriter(out)) {
						print(w,gson,L1);
						}
					}
				}
			else
				{
				final List<JsonArray> L2 = readLog(Paths.get(args.get(1)));
				final List<JsonArray> L3 = new ArrayList<>();
				
				while(!L1.isEmpty() && !L2.isEmpty()) {
					final JsonArray oa = L1.get(0);
					final JsonArray ob = L2.get(0);
					if(equals(oa,ob)) {
						L1.remove(0);
						L2.remove(0);
						L3.add(oa);
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
						L3.add(oa);
						}
					else
						{
						x++;
						}
					}
				
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outFile)) {
					try(JsonWriter w=gson.newJsonWriter(out)) {
						w.beginArray();
						
						
						if(!hide_file1) {
							w.beginObject();
							w.name("name");
							w.value("L1");
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
							w.name("name");
							w.value("L2");
							w.name("description");
							w.value("Unique to "+args.get(1));
							w.name("file");
							w.value(args.get(1));
							w.name("items");
							print(w,gson,L2);
							w.endObject();
							}
						if(!hide_file3) {
							w.beginObject();
							w.name("name");
							w.value("L3");
							w.name("description");
							w.value("common to "+args.get(0)+" and "+args.get(1));
							w.name("file");
							w.beginArray();
							w.value(args.get(0));
							w.value(args.get(1));
							w.endArray();
							w.name("items");
							print(w,gson,L3);
							w.endObject();
							}
						
						if(!hide_unique && L1.size()==1 && L2.size()==1) {
							final JsonArray oa = L1.get(0);
							final JsonArray ob = L2.get(0);
							final Set<JsonObject> hash_1 = new HashSet<>(oa.size());
							for(int i=0;i< oa.size();i++) {
								hash_1.add(oa.get(i).getAsJsonObject());
								}
							final Set<JsonObject> hash_2 = new HashSet<>(ob.size());
							for(int i=0;i< ob.size();i++) {
								hash_2.add(ob.get(i).getAsJsonObject());
								}
							final HashSet<JsonObject> hash_comm = new HashSet<>(hash_1);
							hash_comm.retainAll(hash_2);
							hash_1.removeAll(hash_comm);
							hash_2.removeAll(hash_comm);
							
							w.beginObject();
							w.name("name");
							w.value("uniq");
							w.name("description");
							w.value("common/discordant values shared between the unique item of L1 and the unique item of L2");
							
							w.name("common");
							w.beginArray();
							for(JsonObject o: hash_comm) {
								w.jsonValue(gson.toJson(o));
								}
							w.endArray();
							
							w.name("uniq_1");
							w.beginArray();
							for(JsonObject o: hash_1) {
								w.jsonValue(gson.toJson(o));
								}
							w.endArray();
							
							w.name("uniq_2");
							w.beginArray();
							for(JsonObject o: hash_2) {
								w.jsonValue(gson.toJson(o));
								}
							w.endArray();
							
							
							w.endObject();
							}
						
						w.endArray();
						w.flush();
					}
				}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		}
	
	public static void main(String[] args) {
		new NFHashTools().instanceMainWithExit(args);
	}


}
