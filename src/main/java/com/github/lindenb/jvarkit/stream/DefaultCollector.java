/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.stream;

import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;

/** basic implementation of a Collector */
public class DefaultCollector<T, A, R> implements Collector<T, A, R> {

    private final Supplier<A> supplier;
    private final BiConsumer<A, T> accumulator;
    private final BinaryOperator<A> combiner;
    private final Function<A, R> finisher;
    private final Set<Characteristics> characteristics;
    
	public DefaultCollector(
			final Supplier<A> supplier, 
			final BiConsumer<A, T> accumulator, 
			final BinaryOperator<A> combiner,
			final Function<A, R> finisher, 
			final Set<Characteristics> characteristics
			) {
		this.supplier = supplier;
		this.accumulator = accumulator;
		this.combiner = combiner;
		this.finisher = finisher;
		this.characteristics = characteristics;
		}
	@SuppressWarnings("unchecked")
	public DefaultCollector(
			final Supplier<A> supplier, 
			final BiConsumer<A, T> accumulator, 
			final BinaryOperator<A> combiner,
			final Set<Characteristics> characteristics
			) {
		this(supplier,accumulator,combiner,
			(a)->((R)a),
			characteristics
			);
		}

	@Override
	public Supplier<A> supplier() {
		return supplier;
		}
	@Override
	public BiConsumer<A, T> accumulator() {
		return accumulator;
		}
	@Override
	public BinaryOperator<A> combiner() {
		return combiner;
		}
	@Override
	public Function<A, R> finisher() {
		return this.finisher;
		}
	@Override
	public Set<Characteristics> characteristics() {
		return this.characteristics;
		}

	}
