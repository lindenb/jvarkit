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
package com.github.lindenb.jvarkit.swing;

import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;

public class PropertyChangeObserver<T> {
private T value;
private final PropertyChangeSupport pcs;
public PropertyChangeObserver() {
	this(null);
	}
public PropertyChangeObserver(final T value) {
	this.value = value;
	this.pcs = new PropertyChangeSupport(this);
	}
public boolean isPresent() {
	return this.value!=null;
	}
public void reset() {
	this.setValue(null);
	}
public T getValue() {
	return this.value;
	}
public T orElse(final T t) {
	return this.value==null?t:this.value;
	}
public T orThrow() {
	if(this.value==null) throw new NullPointerException("The is no data in this "+this.getClass());
	return this.value;
	}
protected String getEventName() {
	return "value";
	}
public synchronized T setValue(final T newValue) {
	final T oldValue = this.value;
	this.value = newValue;
	this.pcs.firePropertyChange(getEventName(), oldValue, newValue);
	return oldValue;
	}
public void addPropertyChangeListener(PropertyChangeListener listener) {
    this.pcs.addPropertyChangeListener(listener);
	}

public void removePropertyChangeListener(PropertyChangeListener listener) {
    this.pcs.removePropertyChangeListener(listener);
	}
@Override
public int hashCode() {
	return this.value==null?-1:this.value.hashCode();
	}
@Override
public String toString() {
	return String.valueOf(this.value);
	}
}
