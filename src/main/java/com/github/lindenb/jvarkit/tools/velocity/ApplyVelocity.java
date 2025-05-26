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
package com.github.lindenb.jvarkit.tools.velocity;

import java.io.File;
import java.io.PrintWriter;
import java.io.Reader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.SampleSheet;
import com.github.lindenb.jvarkit.io.SampleSheetFactory;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.velocity.VelocityTools;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFReader;

import org.apache.velocity.Template;
import org.apache.velocity.VelocityContext;
import org.apache.velocity.app.VelocityEngine;
import org.apache.velocity.context.Context;
import org.apache.velocity.runtime.RuntimeConstants;
import org.apache.velocity.runtime.resource.loader.StringResourceLoader;

/**

BEGIN_DOC


text utility using apache velocity templates to generate text.

## Syntax

```
cat template.vm | java -jar jvarkit.jar applyvelocity -m manifest.tsv > out.txt
cat template.vm | java -jar jvarkit.jar applyvelocity -J mycontext file.json > out.txt
java -jar jvarkit.jar applyvelocity -m manifest.tsv template1.vm template2.vm > out.txt
```

## Manifest

TSV file with 3 columns

- key :the name of the context injected in velocity
- type : type of value: 'int', 'long', 'float', 'double', 'string', 'boolean' , 'json' (might be a json string (value starts with '{' or '[' ) or a path to a json file ) , 'dict' htsjdk dictionary,
  'vcf' the variants in a vcf file. 'class' a java class, 'instance-of' instance of given java class
- value: the value 


##Example

```
$ cat jeter.json
{
"message":"world"
}

$ echo 'hello ${j.message}' | java -jar dist/jvarkit.jar applyvelocity -J j jeter.json
hello world
```

END_DOC

 */
@Program(name="applyvelocity",
keywords={"velocity","json"},
description="Execute apache velocity macros ",
creationDate="20241023",
modificationDate="20241023",
jvarkit_amalgamion = true
)
public class ApplyVelocity extends Launcher{
	private static final Logger LOG = Logger.of(ApplyVelocity.class);
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	Path outputFile = null;
    @Parameter(names={"-m","--manifest"},description="TSV file with the following required header: key/value/type")
    private Path manifestPath=null;
    @Parameter(names={"--no-flatten"},description="don't convert google json to standard java types.",hidden = true)
    private boolean no_flatten=false;
    @Parameter(names={"-J","--json-file"},description="Load JSON file. Syntax -J KEY file.json.",arity = 2)
    private List<String> kvJsonFile = new ArrayList<>();
    @Parameter(names={"-V","--vcf-file"},description="Load VCF file. Syntax -V KEY file.vcf.",arity = 2)
    private List<String> kvVcfFile = new ArrayList<>();
    @Parameter(names={"-R","--dict-file"},description="Load DICT file. Syntax -R KEY file.dict.",arity = 2)
    private List<String> kvDictFile = new ArrayList<>();


	private Object flatten(final JsonElement e) {
		if(no_flatten) return e;
		if(e.isJsonNull()) {
			return null;
		} else if(e.isJsonPrimitive()) {
			final JsonPrimitive prim = e.getAsJsonPrimitive();
			if(prim.isBoolean()) {
				return prim.getAsBoolean();
				}
			else if(prim.isNumber()) {
				return prim.getAsNumber();
				}
			else if(prim.isString()) {
				return prim.getAsString();
				}
			else
				{
				throw new IllegalStateException(String.valueOf(prim));
				}
			}
		else if(e.isJsonObject()) {
			final JsonObject x = e.getAsJsonObject();
			final Map<String,Object> w = new LinkedHashMap<>();
			for(Map.Entry<String, JsonElement> KV:x.entrySet()) {
				w.put(KV.getKey(),flatten(x.get(KV.getKey())));
				};
			return w;
		}else if(e.isJsonArray()) {
			final JsonArray x = e.getAsJsonArray();
			final List<Object> w = new ArrayList<>(x.size());
			for(int i=0;i< x.size();++i) {
				w.add(flatten(x.get(i)));
				}
			return w;
		}
		throw new IllegalArgumentException(String.valueOf(e));
	}

	private void put(final Context context,String key,Object o) {
		if(StringUtils.isBlank(key)) {
			throw new IllegalArgumentException("empty key ");
			}
		if(context.containsKey(key)) {
			throw new IllegalArgumentException("duplicate key "+key);
			}
		context.put(key,o);
		}
    
    
    @Override
    public int doWork(List<String> args) {
    
    final Context context = new VelocityContext();
	context.put("tool", new VelocityTools());
  	
	try {
			for(int i=0;i+1< kvJsonFile.size();i+=2) {				
				try(Reader r = Files.newBufferedReader(Paths.get(kvJsonFile.get(i+1)))) {					
					final JsonParser jr = new JsonParser();
					final Object o =  flatten(jr.parse(r));
					System.err.println("hello3:"+o);

					put(context, kvJsonFile.get(i+0),o);
	  				}
			}
		

			for(int i=0;i+1< kvVcfFile.size();i+=2) {				
				try(VCFReader r= new VCFFileReader(Paths.get(kvVcfFile.get(i+1)),false)) {
					try(CloseableIterator<VariantContext> iter=r.iterator()) {
						put(context, kvVcfFile.get(i+0), iter.stream().collect(Collectors.toList()));
						}
					}
			}
		

			for(int i=0;i+1< kvDictFile.size();i+=2) {				
				put(context,
						kvDictFile.get(i+0),
						new SequenceDictionaryExtractor().extractDictionary(Paths.get(kvVcfFile.get(i+1)))
						);
			}
		
		
		if(this.manifestPath!=null) {
			final SampleSheet samplesheet = new SampleSheetFactory().
					splitter(CharSplitter.TAB).
					of(this.manifestPath);
			
			final FileHeader header = samplesheet.getHeader();
			header.assertColumnExists("name");
			header.assertColumnExists("key");
			header.assertColumnExists("value");
			for(final FileHeader.RowMap row: samplesheet) {
				final String key = row.get("key");
				final String type= row.get("type").toLowerCase();
				final String value= row.get("value");
				final Object o;
				if(type.equals("string")) {
					o=value;
					}
				else if(type.equals("bool") || type.equals("boolean")) {
					o = Boolean.parseBoolean(value);
					}
				else if(type.equals("int") || type.equals("integer")) {
					o = Integer.parseInt(value);
					}
				else if(type.equals("long")) {
					o = Long.parseLong(value);
					}
				else if(type.equals("float")) {
					o = Float.parseFloat(value);
					}
				else if(type.equals("double")) {
					o = Double.parseDouble(value);
					}
				else if(type.equals("yaml")) {
					LOG.error("yaml is not supported");
					return -1;
					}
				else if(type.equals("vcf")) {
					try(VCFReader r= new VCFFileReader(Paths.get(value),false)) {
						try(CloseableIterator<VariantContext> iter=r.iterator()) {
							o =  iter.stream().collect(Collectors.toList());
							}
						}
					}
				else if(type.equals("dict")) {
					o  = new SequenceDictionaryExtractor().extractDictionary(Paths.get(value));
					}
				else if(type.equals("json")) {
					final Path f = ((value.trim().startsWith("[") || value.trim().startsWith("{"))?null:Paths.get(value));
					final JsonParser jr = new JsonParser();
					if(f==null || !Files.exists(f)) {
							o=flatten(jr.parse(value));
							}
						else
							{
							try(Reader r = Files.newBufferedReader(f)) {
				  				o = flatten(jr.parse(r));
				  				}
							}
					}
				else if(type.equals("class")) {
					final Class<?> c=Class.forName(value);
						o=c;
					}
				else if(type.equals("instance-of") || type.equals("instanceof")) {
					Class<?> c=Class.forName(value);
						o=c.getConstructor().newInstance();
					}
				else	{
					LOG.error("unknow type:"+type+" in "+row);
					return -1;
					}
				put(context,key,o);
				}
			}
		
		for(final String k: context.getKeys()) {
			System.err.println(""+k+"="+context.get(k));
		}
    try(final PrintWriter w= super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
    	if(args.isEmpty())
    		{
            final VelocityEngine ve = new VelocityEngine();
            ve.setProperty(RuntimeConstants.RESOURCE_LOADER, "string");
            ve.setProperty("string.resource.loader.class","org.apache.velocity.runtime.resource.loader.StringResourceLoader");
            ve.init();
            final String templText = IOUtils.slurpInputStream(System.in);
            StringResourceLoader.getRepository().putStringResource("t1", templText);
            final Template template = ve.getTemplate("t1");
            template.merge( context, w);
    		}
    	else
	    	{
    		for(final String filename: args)
                    {
                    final File file=new File(filename);
                    final VelocityEngine ve = new VelocityEngine();
                    ve.setProperty(RuntimeConstants.RESOURCE_LOADER, "file");
                    ve.setProperty("file.resource.loader.class","org.apache.velocity.runtime.resource.loader.FileResourceLoader");
                    //ve.setProperty("file.resource.loader.path",Collections.singletonList(file.getParent()==null?new File("."):file.getParent()));
                    ve.setProperty("file.resource.loader.path",file.getParent()==null?".":file.getParent().toString());
                    ve.init();
                    final Template template = ve.getTemplate(file.getName());
                    template.merge( context, w);
                    }
	    	}
    	w.flush();
    	}
	return 0;
	}
catch(final Throwable err) {
	LOG.error(err);
	return -1;
	}
}

	public static void main(String[] args) {
		new ApplyVelocity().instanceMainWithExit(args);
	}

}
