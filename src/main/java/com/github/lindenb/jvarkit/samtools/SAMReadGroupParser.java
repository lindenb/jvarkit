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

*/
package com.github.lindenb.jvarkit.samtools;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Function;

import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.SAMReadGroupRecord;

/** a simple class parsing @RG provided by user, when generating SAM */
public class SAMReadGroupParser implements Function<String, SAMReadGroupRecord> {
/** parses a valid read group, throw IllegalArgumentException pn error */
@Override
public SAMReadGroupRecord apply(final String readGroupStr) {
	if(readGroupStr==null) throw new IllegalArgumentException("null read group");
	final String tokens[]=CharSplitter.TAB.split(readGroupStr);
	if(!tokens[0].equals("@RG")) {
		throw new IllegalArgumentException("Doesn't start with a @RG / bad read group "+readGroupStr);
		}
	final Map<String,String> props= new LinkedHashMap<>(tokens.length);
	for(int i=0;i< tokens.length;i++)
		{
		final int colon = tokens[i].indexOf(':');
		if(colon<=0) {
			throw new IllegalArgumentException("Colon missing in "+tokens[i]);
			}
		final String key = tokens[i].substring(colon+1);
		if(props.containsKey(key)) throw new IllegalArgumentException("duplicate key "+key+" in "+readGroupStr);
		props.put(key, tokens[i].substring(colon+1));
		}
	final String readGroupId = props.get("ID");
	if(StringUtils.isBlank(readGroupId)) {
		throw new IllegalArgumentException("no ID: in  "+readGroupStr);
		}
	props.remove("ID");
	final SAMReadGroupRecord rg = new SAMReadGroupRecord(readGroupId);
	props.entrySet().forEach(KV->{
		rg.setAttribute(KV.getKey(), KV.getValue());
		});
	return rg;
	}
}
