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
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Supplier;


public class AutoMap<K,V,C> extends AbstractMap<K,C> {
	private BiFunction<K,V,C> collectionMaker;
	private final Map<K,C> delegate;
	private final BiConsumer<C,V> inserter;
	public AutoMap(
			Supplier<Map<K,C>> mapMaker,
			BiFunction<K,V,C> collectionMaker,
			BiConsumer<C,V> inserter
			) {
		this.delegate = mapMaker.get();
		this.collectionMaker = collectionMaker;
		this.inserter=inserter;
		}
	
	public AutoMap(
			BiFunction<K,V,C> collectionMaker,
			BiConsumer<C,V> inserter
			) {
		this(()->new HashMap<K,C>(),collectionMaker,inserter);
		}
	
	public C insert(K key) {
		C col = this.delegate.get(key);
		if(col==null) {
			col = Objects.requireNonNull(this.collectionMaker.apply(key,null));
			this.delegate.put(key, col);
			}
		return col;
		}

	
	public C insert(K key, V value) {
		Objects.requireNonNull(value,"value is null");
		C col = this.delegate.get(key);
		if(col==null) {
			col = Objects.requireNonNull(this.collectionMaker.apply(key,value));
			this.delegate.put(key, col);
			}
		this.inserter.accept(col, value);
		return col;
		}
	
	public static <K,V> AutoMap<K,V,List<V>> makeList() {
		return new AutoMap<K,V,List<V>>(
				()->new HashMap<>(),
				(A,B)->new ArrayList<V>(),
				(A,B)->A.add(B)
				);
		}
	public static <K,V> AutoMap<K,V,Set<V>> makeSet() {
		return new AutoMap<K,V,Set<V>>(
				()->new HashMap<>(),
				(A,B)->new HashSet<V>(),
				(A,B)->A.add(B)
				);
		}

	/** create a simple key value */
	public static <K,V> AutoMap<K,V,V> make( BiFunction<K,V,V> collectionMaker) {
		return new AutoMap<K,V,V>(
				()->new HashMap<>(),
				(A,B)->collectionMaker.apply(A,B),
				(A,B)->{}
				);
		}
	/** create a simple key value , with a supplier that doesn't need the key to be instanied */
	public static <K,V> AutoMap<K,V,V> make( Supplier<V> simpleValueMaker) {
		return new AutoMap<K,V,V>(
				()->new HashMap<>(),
				(A,B)->simpleValueMaker.get(),
				(A,B)->{}
				);
		}
	
	public int size() {
		return delegate.size();
	}

	public boolean isEmpty() {
		return delegate.isEmpty();
	}

	public boolean containsKey(Object key) {
		return delegate.containsKey(key);
	}

	public boolean containsValue(Object value) {
		return delegate.containsValue(value);
	}

	public C get(Object key) {
		return delegate.get(key);
	}

	public C put(K key, C value) {
		return delegate.put(key, value);
	}

	public C remove(Object key) {
		return delegate.remove(key);
	}

	public void putAll(Map<? extends K, ? extends C> m) {
		delegate.putAll(m);
	}

	public void clear() {
		delegate.clear();
	}

	public Set<K> keySet() {
		return delegate.keySet();
	}

	public Collection<C> values() {
		return delegate.values();
	}

	public Set<Entry<K, C>> entrySet() {
		return delegate.entrySet();
		}



	
	}
