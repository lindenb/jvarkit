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

import java.util.Collection;
import java.util.Objects;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import com.sleepycat.bind.tuple.TupleBinding;
import com.sleepycat.bind.tuple.TupleInput;
import com.sleepycat.bind.tuple.TupleOutput;

/** Berkeleydb binding for java.util.collection */
public class CollectionBinding<T extends Collection<TYPE>,TYPE> extends TupleBinding<T> {
	private final IntFunction<T> generator;
	private final TupleBinding<TYPE> binding;
	
	/**
	 * 
	 * @param generator how to create a collection with a reserve of 'x' items
	 * @param binding how to serialize TYPE
	 */
	public CollectionBinding(final IntFunction<T> generator,TupleBinding<TYPE> binding) {
		this.generator= Objects.requireNonNull(generator);
		this.binding = Objects.requireNonNull(binding);
		}
	/**
	 * 
	 * @param generator how to create a collection
	 * @param binding how to serialize TYPE
	 */
	public CollectionBinding(final Supplier<T> generator,TupleBinding<TYPE> binding) {
		this(I->generator.get(),binding);
		}
	
	@Override
	public T entryToObject(final TupleInput in) {
		final int n =  in.readInt();
		final T col = this.generator.apply(n);
		for(int i=0;i< n;++i) {
			col.add(binding.entryToObject(in));
			}
		return col;
		}

	@Override
	public void objectToEntry(T t, TupleOutput out) {
		out.writeInt(t.size());
		for(TYPE item: t) {
			binding.objectToEntry(item,out);
			}
		}
		
	
}
