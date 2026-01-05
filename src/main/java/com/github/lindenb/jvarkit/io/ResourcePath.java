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
package com.github.lindenb.jvarkit.io;

import java.io.IOException;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;

import com.beust.jcommander.IStringConverter;
import com.github.lindenb.jvarkit.lang.AttributeMap;
import com.github.lindenb.jvarkit.lang.StringUtils;

/** paths associated to attributes */
public interface ResourcePath extends AttributeMap {
public static final String OPT_DESC="A file followed with attributes. Syntax: `path/to/file(;key[:=]value)*`  . ";

public static class StringConverter implements IStringConverter<ResourcePath> {
	@Override
	public ResourcePath convert(final String s) {
		return ResourcePath.of(s);
		}
	}

/** get Path */
public Path getPath();

/** does file exist */
public default boolean exists() {
	return Files.exists(getPath());
}

/** implementation of ResourcePath */
static class ResourcePathImpl implements ResourcePath {
	String declaration;
	Path path;
	final Map<String, String> properties = new LinkedHashMap<>();
	@Override
	public Path getPath() {
		return this.path;
		}
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
	rsrc.path = Paths.get(s.substring(0,x1).trim());
	s=s.substring(x1+1);
	if(!StringUtils.isBlank(s)) {
		try(StringReader bais = new StringReader(s)) {
			final StreamTokenizer st = new StreamTokenizer(bais) {
				@Override
				public void parseNumbers() {
					// https://stackoverflow.com/questions/46536511/
					}
				};
			st.wordChars('_', '_');
			st.wordChars('-','-');
			st.wordChars('0', '9');
			st.slashSlashComments(false);
			st.slashStarComments(false);
			st.eolIsSignificant(false);
			for(;;) {
				int op =st.nextToken();
				if(op==StreamTokenizer.TT_EOF) break;
				final String key;
				if(op==StreamTokenizer.TT_NUMBER) throw new IllegalArgumentException();
				if(op==';') continue;//empy decl
				if(op=='"' || op=='\'' || op==StreamTokenizer.TT_WORD) {
					key = st.sval;
					}
				else
					{
					throw new IllegalArgumentException("key expected");
					}
				if(rsrc.properties.containsKey(key)) {
					throw new IllegalArgumentException("duplicate key "+key+" in "+rsrc.declaration);
					}
				op =st.nextToken();
				if(op==StreamTokenizer.TT_EOF || op==';') {
					rsrc.properties.put(key,"true");
					if(op==StreamTokenizer.TT_EOF) break;
					continue;
					}
				if(!(op==':' || op=='=')) {
					throw new IllegalArgumentException("expected ':' or '=' after key="+key+" in "+rsrc.declaration);
					}
				final String value;
				op =st.nextToken();
				if(op==StreamTokenizer.TT_NUMBER) throw new IllegalArgumentException();
				if(op=='"' || op=='\'' || op==StreamTokenizer.TT_WORD) {
					value = st.sval;
					}
				else {
					throw new IllegalArgumentException("value expected for key="+key+" in "+rsrc.declaration);
					}
				rsrc.properties.put(key,value);
				op =st.nextToken();
				if(op==StreamTokenizer.TT_EOF) break;
				if(op!=';') {
					throw new IllegalArgumentException("expected ';' key="+key+"="+value+" in "+rsrc.declaration);
					}
				}//end for
			}//end string reader
		catch(final IOException err) {
			throw new IllegalArgumentException(err);
			}
		}//end if 
	return rsrc;
	}
}
