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
package com.github.lindenb.jvarkit.util;

import java.util.AbstractMap;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class FunctionalMap<K,V> implements Iterable<Map.Entry<K, V>>{
	private final Map<K,V> delegate;
	public FunctionalMap() {
		delegate = new HashMap<>();
		}
	public FunctionalMap(final FunctionalMap<K,V> cp) {
		this(cp.delegate);
		}
	public FunctionalMap(final Map<K,V> hash) {
		delegate = new HashMap<>(hash);
		}
	
	public FunctionalMap<K,V> plus(K k1,V v1) {
		final FunctionalMap<K,V> cp = makeCopy();
		cp.delegate.put(k1, v1);
		return cp;
		}
	
	
	public FunctionalMap<K,V> minus(K k1) {
		if(!containsKey(k1)) return this;
		final FunctionalMap<K,V> cp = makeCopy();
		cp.delegate.remove(k1);
		return cp;
		}
	
	@Override
	public int hashCode() {
		return delegate.hashCode();
		}
	@Override
	public boolean equals(Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof FunctionalMap)) return false;
		return delegate.equals( FunctionalMap.class.cast(obj).delegate);
		}
	
	public Set<K> keySet() {
		return Collections.unmodifiableSet(this.delegate.keySet());
		}
	public Collection<V> values() {
		return Collections.unmodifiableCollection(this.delegate.values());
		}
	
	public Set<Map.Entry<K, V>> entrySet() {
		return stream().collect(Collectors.toSet());
		}
	
	public Stream<Map.Entry<K, V>> stream() {
		return this.delegate.entrySet().
				stream().
				/** be sure that delegate is inmutable */
				map(KV->new AbstractMap.SimpleEntry<K,V>(KV.getKey(),KV.getValue()))
				;
		}
	
	@Override
	public Iterator<Entry<K, V>> iterator() {
		return stream().iterator();
		}
	
	public V get(K k) {
		return this.delegate.get(k);
		}
	
	public V getOrDefaul(K k,V def) {
		return this.delegate.getOrDefault(k, def);
		}
	
	public boolean containsKey(final K k) {
		return this.delegate.containsKey(k);
		}
	
	@Override
	public String toString() {
		return delegate.toString();
		}
	
	public int size() {
		return this.delegate.size();
		}
	
	protected FunctionalMap<K,V> makeCopy() {
		return new FunctionalMap<>(this);
		}
	@Override
	public FunctionalMap<K,V> clone() {
		return this;
		}
	}
