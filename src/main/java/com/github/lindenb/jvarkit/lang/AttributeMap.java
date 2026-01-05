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
package com.github.lindenb.jvarkit.lang;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.OptionalInt;
import java.util.OptionalLong;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/** general interface to implement a Map<String,String> */
public interface AttributeMap {

/** get all attributes */
public Map<String,String> getAttributes();
/** get attribute as string or return default */
public default String getAttribute(final String key,String defaultValue) {
	return getAttribute(key).orElse(defaultValue);
	}
/** get Optional attribute */
public default Optional<String> getAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) {
		reportMissingKey(key);
		return Optional.empty();
		}
	return Optional.of(m.get(key));
	}
/** return true of map has attribute */
public default boolean hasAttribute(final String key) {
	return getAttributes().containsKey(key);
	}


/** get Optional attribute of enum item , ignore case */
public default  <E extends Enum<E>>  Optional<E> getEnumAttribute(final Class<E> theEnum,final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) {
		reportMissingKey(key);
		return Optional.empty();
		}
	final String s = m.get(key).trim();
	for (E c : java.util.EnumSet.allOf(theEnum)) {
		 if(c.name().equalsIgnoreCase(s)) return Optional.of(c);
		 }
	throw new IllegalArgumentException("Cannot convert "+key+"=\""+s+"\" to items: "+ 
			java.util.EnumSet.allOf(theEnum).stream().
				map(S->S.name()).
				collect(Collectors.joining(" | ")));
	}

/** get Optional int attribute */
public default OptionalInt getIntAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) {
		reportMissingKey(key);
		return OptionalInt.empty();
		}
	try {
		return OptionalInt.of(Integer.parseInt(m.get(key)));
		} catch(final NumberFormatException err) {
			throw new IllegalStateException("Cannot parse int attribute for key \""+key+"\".",err);
		}
	}

/** get Optional long attribute */
public default OptionalLong getLongAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) {
		reportMissingKey(key);
		return OptionalLong.empty();
		}
	try {
		return OptionalLong.of(Long.parseLong(m.get(key)));
		} catch(final NumberFormatException err) {
			throw new IllegalStateException(key,err);
		}
	}

/** get Optional double attribute */
public default OptionalDouble getDoubleAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) {
		reportMissingKey(key);
		return OptionalDouble.empty();
		}
	try {
		return OptionalDouble.of(Double.parseDouble(m.get(key)));
		} catch(final NumberFormatException err) {
			throw new IllegalStateException(key,err);
		}
	}

public default boolean getBooleanAttribute(final String key) {
	final Map<String,String> m = getAttributes();
	if(!m.containsKey(key)) {
		reportMissingKey(key);
		return false;
		}
	final String val = m.get(key);
	if(val.equalsIgnoreCase("false")) return false;
	if(val.equalsIgnoreCase("no")) return false;
	if(val.equalsIgnoreCase("true")) return true;
	if(val.equalsIgnoreCase("yes")) return true;
	throw new IllegalStateException("not a boolean "+key+":"+val);
	}


/** get all pairs(key/value) */
public default Stream<Map.Entry<String, String>> entries() {
	return getAttributes().entrySet().stream();
	}

/** give a chance to report a missing key */
public default void reportMissingKey(final String keyName) {
	}

static String _toString(final Map<String,String> map)  {
	final StringBuilder sb=new StringBuilder();
	int i=1;
	for(final String key: map.keySet()) {
		sb.append("$").append(i++).append(" : ").
			append(key).append(" : ").
			append(map.get(key)).append("\n");
		}
	return sb.toString();
	}

/** wraps a java.util.map */
public static AttributeMap wrap(final Map<String, String> hash) {
	final Map<String,String> umap = Collections.unmodifiableMap(hash);
	return new AttributeMap() {
		@Override
		public Map<String, String> getAttributes() {
			return umap;
			}
		@Override
		public String toString() {
			return _toString(this.getAttributes());
			}
		};
	}

/** wraps a AttributeMap reporting a missing key */
public static AttributeMap verbose(final AttributeMap src,final Consumer<String> reportMissingKey) {
	if(reportMissingKey==null) return src;
	return new AttributeMap() {
		final Set<String> seen = new HashSet<>();
		@Override
		public Map<String, String> getAttributes() {
			return src.getAttributes();
			}
		@Override
		public void reportMissingKey(String keyName) {
			if(seen.add(keyName)) {
				reportMissingKey.accept(keyName);
				}
			}
		@Override
		public String toString() {
			return _toString(this.getAttributes());
			}
		};
	}
/** create a non-null empty attributeMap */
public static AttributeMap empty() {
	return wrap(Collections.emptyMap());
	}
/** create from a pairs of string k1,v1,k2,v2,k3,v3... */
public static AttributeMap fromPairs(String...array) {
	if(array.length==0) return empty();
	if(array.length%2!=0) throw new IllegalArgumentException("not an odd number of strings");
	final Map<String,String> hash = new HashMap<>(array.length/2);
	for(int i=0;i+1< array.length;i+=2) {
		hash.put(array[i], array[i+1]);
		}
	return wrap(hash);
	}
}
