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
