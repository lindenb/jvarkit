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
package com.github.lindenb.jvarkit.bdb;

import java.util.Map;
import java.util.Objects;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

/** Berkeleydb binding for java.util.Map */
public class MapBinding<T extends Map<KEY,VALUE>,KEY,VALUE> extends TupleBinding<T> {
	private final IntFunction<T> generator;
	private final TupleBinding<KEY> keyBinding;
	private final TupleBinding<VALUE> valueBinding;
	
	/**
	 * 
	 * @param generator how to create a Map with a reserve of 'x' items
	 * @param keyBinding how to serialize the key
	 * @param valueBinding how to serialize the value
	 */
	public MapBinding(final IntFunction<T> generator,
			final TupleBinding<KEY> keyBinding,
			final TupleBinding<VALUE> valueBinding
			) {
		this.generator= Objects.requireNonNull(generator);
		this.keyBinding = Objects.requireNonNull(keyBinding);
		this.valueBinding = Objects.requireNonNull(valueBinding);
		}
	/**
	 * 
	 * @param generator how to create a collection
	 * @param keyBinding how to serialize the key
	 * @param valueBinding how to serialize the value
	 */
	public MapBinding(final Supplier<T> generator,			
			final TupleBinding<KEY> keyBinding,
			final TupleBinding<VALUE> valueBinding) {
		this(I->generator.get(),keyBinding,valueBinding);
		}
	
	@Override
	public T entryToObject(final TupleInput in) {
		final int n =  in.readInt();
		final T col = this.generator.apply(n);
		for(int i=0;i< n;++i) {
			final KEY k = this.keyBinding.entryToObject(in);
			final VALUE v = this.valueBinding.entryToObject(in);
			col.put(k, v);
			}
		return col;
		}

	@Override
	public void objectToEntry(final T t, TupleOutput out) {
		out.writeInt(t.size());
		for(KEY key: t.keySet()) {
			this.keyBinding.objectToEntry(key, out);
			this.valueBinding.objectToEntry(t.get(key), out);
			}
		}
		
	
}
