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
package com.github.lindenb.jvarkit.lang;

import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.OptionalLong;
import java.util.stream.Stream;

/** general interface to implement a Map<String,String> */
public interface AttributeMap {

/** get all attributes */
public Map<String,String> getAttributes();
/** get attribute as string or return default */
public default String getAttribute(final String key,String defaultValue) {
	return getAttributes().getOrDefault(key,defaultValue);
	}
/** get Optional attribute */
public default Optional<String> getAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) return Optional.empty();
	return Optional.of(m.get(key));
	}
/** return true of map has attribute */
public default boolean hasAttribute(final String key) {
	return getAttributes().containsKey(key);
	}

/** get Optional int attribute */
public default OptionalInt getIntAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) return OptionalInt.empty();
	try {
		return OptionalInt.of(Integer.parseInt(m.get(key)));
		} catch(final NumberFormatException err) {
			throw new IllegalStateException(key,err);
		}
	}

/** get Optional long attribute */
public default OptionalLong getLongAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) return OptionalLong.empty();
	try {
		return OptionalLong.of(Long.parseLong(m.get(key)));
		} catch(final NumberFormatException err) {
			throw new IllegalStateException(key,err);
		}
	}

/** get Optional double attribute */
public default OptionalDouble getDoubleAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) return OptionalDouble.empty();
	try {
		return OptionalDouble.of(Double.parseDouble(m.get(key)));
		} catch(final NumberFormatException err) {
			throw new IllegalStateException(key,err);
		}
	}

public default boolean getBooleanAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) return false;
	final String val = m.get(key);
	if(val.equalsIgnoreCase("false")) return false;
	if(val.equalsIgnoreCase("true")) return true;
	throw new IllegalStateException("not a boolean "+key+":"+val);
	}


/** get all pairs(key/value) */
public default Stream<Map.Entry<String, String>> entries() {
	return getAttributes().entrySet().stream();
	}
}
