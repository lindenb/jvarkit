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

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Optional;

import com.github.lindenb.jvarkit.lang.StringUtils;

public interface ResourcePath {

public Map<String,String> getAttributes();
		
public default String getAttribute(final String key,String defaultValue) {
	return getAttribute(key).orElse(defaultValue);
	}

public default Optional<String> getAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) return Optional.empty();
	return Optional.of(m.get(key));
	}

public default boolean hasAttribute(final String key) {
	return getAttributes().containsKey(key);
	}

	
static class ResourcePathImpl implements ResourcePath {
	String declaration;
	Path path;
	final Map<String, String> properties = new LinkedHashMap<>();
	@Override
	public Map<String, String> getAttributes() {
		return Collections.unmodifiableMap(properties);
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
	while(!s.isEmpty()) {
		int i=0;
		while(i<s.length() && s.charAt(i)!='=' && s.charAt(i)!=';') {
			i++;
			}
		String key = s.substring(0,i).trim();
		}
	return rsrc;
	}
}
