/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.io;

import java.io.StreamTokenizer;
import java.io.StringReader;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Supplier;

import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;

/** paths associated to attributes */
public interface ResourcePath extends AttributeMap {

/** implementation of ResourcePath */
static class ResourcePathImpl implements ResourcePath {
	String declaration;
	Path path;
	final Map<String, String> properties = new LinkedHashMap<>();
	@Override
	public Map<String, String> getAttributes() {
		return Collections.unmodifiableMap(properties);
		}
	@Override
	public String toString() {
		return declaration;
		}
	}
	
/** create ResourcePath from string */
public static ResourcePath of(String s) {
	if(StringUtils.isBlank(s)) throw new IllegalArgumentException("empty string");
	final ResourcePathImpl rsrc = new ResourcePathImpl();
	rsrc.declaration = s;
	int x1 = s.indexOf(";");
	if(x1==0)  throw new IllegalArgumentException("starts with ';' "+s);
	rsrc.path = Paths.get(s.substring(0,x1));
	s=s.substring(x1+1);
	try(StringReader bais = new StringReader(s)) {
		final StreamTokenizer st = new StreamTokenizer(bais);
		st.wordChars('_', '_');
		st.wordChars('-','-');
		st.wordChars('0', '9');
		st.slashSlashComments(false);
		st.slashStarComments(false);
		final Supplier<String> nextToken= ()-> {
			for(;;)
				final int op=st.nextToken();
				if(op==StreamTokenizer.TT_EOF) return null;
				return st.sval;
				}
			};
		boolean got_colon = false;
		String key = null;
		while((op=st.nextToken())!=StreamTokenizer.TT_EOF) {
			if( op==StreamTokenizer.TT_WORD ||
				op==StreamTokenizer.TT_NUMBER || 
				op == '\'' || op=='\"') {
				String sval ;
				sval = st.sval;
					
				if(key==null) {
					key = sval;
					if(rsrc.properties.containsKey(key)) {
						throw new IllegalArgumentException("duplicate key "+key+" in "+rsrc.declaration);
						}
					got_colon=false;
					}
				else if(got_colon)
					{
					rsrc.properties.put(key,sval);
					key=null;
					got_colon=false;
					}
				else
					{
					throw new IllegalArgumentException("colon missing");
					}
				}
			else if(op=='=' || op==':') {
				got_colon = true;
				}
			else if(op==';') {
				got_colon = false;
				key=null;
				}
			else {
				
				}
			}
		}
	catch(final Throwable err) {
		throw new IllegalArgumentException(err);
		}
	return rsrc;
	}
}
