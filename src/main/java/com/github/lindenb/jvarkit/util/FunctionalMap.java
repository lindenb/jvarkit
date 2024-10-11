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

/**
 * inmmutable Map
 * @author lindenb
 *
 * @param <K> key
 * @param <V> value
 */
public class FunctionalMap<K,V> implements Iterable<Map.Entry<K, V>>{
	private final Map<K,V> delegate;
	public FunctionalMap() {
		delegate = new HashMap<>();
		}
	
	public static <K,V>  FunctionalMap<K,V> make() {
		return new  FunctionalMap<K,V>();
		}
	
	public static <K,V>  FunctionalMap<K,V> of(K k,V v) {
		return new  FunctionalMap<K,V>().plus(k,v);
		}
	public static <K,V>  FunctionalMap<K,V> of(K k1,V v1,K k2,V v2) {
		return of(k1,v1).plus(k2,v2);
		}
	public static <K,V>  FunctionalMap<K,V> of(K k1,V v1,K k2,V v2,K k3,V v3) {
		return of(k1,v1,k2,v2).plus(k3,v3);
		}
	public static <K,V>  FunctionalMap<K,V> of(K k1,V v1,K k2,V v2,K k3,V v3,K k4,V v4) {
		return of(k1,v1,k2,v2,k3,v3).plus(k4,v4);
		}
	public static <K,V>  FunctionalMap<K,V> of(K k1,V v1,K k2,V v2,K k3,V v3,K k4,V v4,K k5,V v5) {
		return of(k1,v1,k2,v2,k3,v3,k4,v4).plus(k5,v5);
		}
	public FunctionalMap(final FunctionalMap<K,V> cp) {
		this(cp.delegate);
		}
	public FunctionalMap(final Map<K,V> hash) {
		delegate = new HashMap<>(hash);
		}
	
	public FunctionalMap(K k1,V v1) {
		this();
		delegate.put(k1, v1);
		}
	public FunctionalMap<K,V> plus(FunctionalMap<K,V> fm) {
		return plus(fm.delegate);
		}
	public FunctionalMap<K,V> plus(final Map<K,V> hash) {
		if(hash.isEmpty()) return this;
		final FunctionalMap<K,V> cp = makeCopy();
		for(final K k:hash.keySet()) {
			cp.delegate.put(k, hash.get(k));
			}
		return cp;
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
	
	public FunctionalMap<K,V> minus(Set<K> set) {
		if(set.stream().noneMatch(key->delegate.containsKey(key))) return this;
		final FunctionalMap<K,V> cp = makeCopy();
		for(final K k:set) {
			cp.delegate.remove(k);
			}
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
	
	public V getOrDefault(K k,V def) {
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
	public boolean isEmpty() {
		return this.delegate.isEmpty();
		}
	
	protected FunctionalMap<K,V> makeCopy() {
		return new FunctionalMap<>(this);
		}
	@Override
	public FunctionalMap<K,V> clone() {
		return this;
		}
	}
