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
package com.github.lindenb.jvarkit.function;

import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.Vector;
import java.util.function.Supplier;
/**
 * A chain of <code>Supplier&lt;Optional&lt;T&gt;&gt;</code>
 */
public class SupplierChain<T> implements Supplier<Optional<T>> {
private final List<Supplier<Optional<T>>> chain = new Vector<>();
public SupplierChain(final Supplier<Optional<T>> one) {
	this.chain.add(Objects.requireNonNull(one));
	}

public SupplierChain(final List<Supplier<Optional<T>>> list) {
	this.chain.addAll(list);
	}
public SupplierChain() {
}
/** add this supplier at the head of the chain  (high priority ) */
public SupplierChain<T> before(final Supplier<Optional<T>> supplier) {
	this.chain.add(0, Objects.requireNonNull(supplier, "SupplierChain::before(null"));
	return this;
	}
/** add this supplier at the tail of the chain  (low priority ) */
public SupplierChain<T> after(final Supplier<Optional<T>> supplier) {
	this.chain.add(Objects.requireNonNull(supplier, "SupplierChain::after(null"));
	return this;
	}

@Override
/** return the 'T' value, never null.  Optional.empty() if no item was found in the chain found something */
public Optional<T> get() {
	return this.chain.stream().
				map(S->S.get()).
				filter(O->Objects.requireNonNull(O,"element in chain returned null").isPresent()).
				findFirst().
				orElse(Optional.empty());
	}

@Override
public String toString() {
	return "SupplierChain<>("+this.chain.size()+")";
	}
}
