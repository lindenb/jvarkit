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
package com.github.lindenb.jvarkit.json;

import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Supplier;

import com.github.lindenb.jvarkit.lang.StringUtils;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;

import htsjdk.samtools.util.IOUtil;

/**
 * utility to convert stream of JsonReader to simple java object (velocity... )
 *
 */
public class FromJson {
	
	private Function<String, Object> numberMapper = S ->{
        if(StringUtils.isInteger(S)) return Integer.parseInt(S);
        if(StringUtils.isLong(S)) return Long.parseLong(S);
        if(StringUtils.isDouble(S)) return Double.parseDouble(S);
        return S;
		};
	private Supplier<List<Object>> listSupplier = () -> new ArrayList<>();
	private Supplier<Map<String,Object>> mapSupplier = () -> new LinkedHashMap<>();
	
	public Object parseJsonPath(final Path path) throws IOException {
		try(JsonReader jsr = new JsonReader(IOUtil.openFileForBufferedReading(path))) {
			return parse(jsr);
			}
		}
	
	public Object parseJsonString(final String str) throws IOException {
		try(JsonReader jsr = new JsonReader(new StringReader(str))) {
			return parse(jsr);
			}
		}
	
	/** set the class generating the list for a json-array */
	public FromJson setListSupplier(Supplier<List<Object>> listSuppliter) {
		this.listSupplier = listSuppliter;
		return this;
		}
	/** set the class generating the map for a json-object */
	public FromJson setMapSupplier(Supplier<Map<String, Object>> mapSupplier) {
		this.mapSupplier = mapSupplier;
		return this;
		}
	
	public Object parse(final JsonReader reader) throws IOException {
        final JsonToken jsonToken = reader.peek();
        switch( jsonToken ) {
            case NULL:
                reader.nextNull();
                return null;
            case STRING:
                return reader.nextString();
            case BOOLEAN:
                return reader.nextBoolean();
            case NUMBER:
                return numberMapper.apply(reader.nextString());
            case BEGIN_ARRAY:
                final List<Object> array = this.listSupplier.get();
                reader.beginArray();
                while( reader.peek() != JsonToken.END_ARRAY )
                    array.add(parse(reader));
                reader.endArray();
                return array;
            case BEGIN_OBJECT:
                final Map<String,Object> object = this.mapSupplier.get();
                reader.beginObject();
                while( reader.peek() != JsonToken.END_OBJECT ) {
                    final String key = reader.nextName();
                    final Object value = parse(reader);
                    object.put(key, value);
                	}
                reader.endObject();
                return object;

            default:
                throw new IllegalStateException("Unexpected JSON token: "+jsonToken);
        	}
		}
}
