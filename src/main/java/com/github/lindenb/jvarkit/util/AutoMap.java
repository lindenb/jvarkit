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
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * map inserting key,value if it doesn't already exists
 * @author lindenb
 *
 * @param <K>
 * @param <V>
 * @param <CONTAINER_OF_V>
 */
public class AutoMap<K,V,CONTAINER_OF_V> extends AbstractMap<K,CONTAINER_OF_V> {
	private BiFunction<K,V,CONTAINER_OF_V> collectionMaker;
	private final Map<K,CONTAINER_OF_V> delegate;
	private final BiConsumer<CONTAINER_OF_V,V> inserter;
	public AutoMap(
			Supplier<Map<K,CONTAINER_OF_V>> mapMaker,
			BiFunction<K,V,CONTAINER_OF_V> collectionMaker,
			BiConsumer<CONTAINER_OF_V,V> inserter
			) {
		this.delegate = mapMaker.get();
		this.collectionMaker = collectionMaker;
		this.inserter=inserter;
		}
	
	public AutoMap(
			BiFunction<K,V,CONTAINER_OF_V> collectionMaker,
			BiConsumer<CONTAINER_OF_V,V> inserter
			) {
		this(()->new HashMap<K,CONTAINER_OF_V>(),collectionMaker,inserter);
		}
	
	public CONTAINER_OF_V insert(K key) {
		CONTAINER_OF_V col = this.delegate.get(key);
		if(col==null) {
			col = Objects.requireNonNull(this.collectionMaker.apply(key,null));
			this.delegate.put(key, col);
			}
		return col;
		}

	
	public CONTAINER_OF_V insert(final K key,final V value) {
		Objects.requireNonNull(value,"value is null");
		CONTAINER_OF_V col = this.delegate.get(key);
		if(col==null) {
			col = Objects.requireNonNull(this.collectionMaker.apply(key,value));
			this.delegate.put(key, col);
			}
		this.inserter.accept(col, value);
		return col;
		}
	
	/** create a new AutoMap<KEY,VALUE,ArrayList<VALUE>> */
	public static <K,V> AutoMap<K,V,List<V>> makeList() {
		return new AutoMap<K,V,List<V>>(
				()->new HashMap<>(),
				(A,B)->new ArrayList<V>(),
				(A,B)->A.add(B)
				);
		}
	/** create a new AutoMap<KEY,VALUE,HashSet<VALUE>> */
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
	
	/** create a simple key value using a function that creates a value from a key*/
	public static <K,V> AutoMap<K,V,V> make( final Function<K,V> createValueFromKey) {
		final BiFunction<K,V,V> collectionMaker = (A,B)->createValueFromKey.apply(A);
		return make(collectionMaker);
		}

	
	/** create a simple key value , with a supplier that doesn't need the key to be instanied */
	public static <K,V> AutoMap<K,V,V> make( Supplier<V> simpleValueMaker) {
		return new AutoMap<K,V,V>(
				()->new HashMap<>(),
				(A,B)->simpleValueMaker.get(),
				(A,B)->{}
				);
		}
	@Override	
	public int size() {
		return delegate.size();
	}

	@Override	
	public boolean isEmpty() {
		return delegate.isEmpty();
	}

	@Override	
	public boolean containsKey(Object key) {
		return delegate.containsKey(key);
	}

	@Override	
	public boolean containsValue(Object value) {
		return delegate.containsValue(value);
	}

	@Override	
	public CONTAINER_OF_V get(Object key) {
		return delegate.get(key);
	}

	@Override	
	public CONTAINER_OF_V put(K key, CONTAINER_OF_V value) {
		return delegate.put(key, value);
	}

	@Override	
	public CONTAINER_OF_V remove(Object key) {
		return delegate.remove(key);
	}

	@Override	
	public void putAll(Map<? extends K, ? extends CONTAINER_OF_V> m) {
		delegate.putAll(m);
	}
	@Override	
	public void clear() {
		delegate.clear();
	}
	@Override	
	public Set<K> keySet() {
		return delegate.keySet();
	}
	@Override	
	public Collection<CONTAINER_OF_V> values() {
		return delegate.values();
	}
	@Override	
	public Set<Entry<K, CONTAINER_OF_V>> entrySet() {
		return delegate.entrySet();
		}

	@Override
	public String toString() {
		return delegate.toString();
		}

	
	}
