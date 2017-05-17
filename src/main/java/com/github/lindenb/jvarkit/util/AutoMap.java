/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.util;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.function.Supplier;

/*** auto filled map */
public class AutoMap<K,V>
	{
	private final Map<K,V> delegate;
	private final Supplier<V> defaultSupplier;

	/** 
	 * 
	 * @param delegate the delegate map
	 * @param defaultSupplier default value supplier for 'get(key)'
	 */
	public AutoMap(final Map<K,V> delegate,final Supplier<V> defaultSupplier) {
		this.delegate = delegate;
		this.defaultSupplier = defaultSupplier;
		}
	/**
	 * auto map backed with a HashMap
	 */
	public AutoMap(final Supplier<V> defaultSupplier)
		{
		this(new HashMap<K,V>(),defaultSupplier);
		}
	
	protected Map<K,V> getDelegate() { return this.delegate;}
	protected Supplier<V> getDefaultValueSupplier() { return this.defaultSupplier;}
	
	public int size() {
		return getDelegate().size();
	}
	public boolean isEmpty() {
		return getDelegate().isEmpty();
	}
	/** put a default value with 'key' */
	public V put(final K key) {
		return this.put(key, this.getDefaultValueSupplier().get());
	}
	
	/** call delegate key, if the key doesn't exist, create a new value with defaultSupplier, insert into
	 * the delegate map, and return the new value
	 * @param key
	 * @return map[key]
	 */
	public V get(final K key) {
		if(getDelegate().containsKey(key))
			{
			return this.getDelegate().get(key);
			}
		else
			{
			V v = this.getDefaultValueSupplier().get();
			this.getDelegate().put(key,v);
			return v;
			}
		}

	
	public V put(final K key,final V value) {
		return getDelegate().put(key, value);
	}
	public Set<K> keySet() {
		return getDelegate().keySet();
	}
	public Collection<V> values() {
		return getDelegate().values();
	}
	
	@Override
	public String toString() {
		return getDelegate().toString();
		}
	}