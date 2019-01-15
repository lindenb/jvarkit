package com.github.lindenb.jvarkit.table;

import java.util.function.Function;

public interface Column {
public String getName();
public void setStringConverter(final Function<Object, String> fun);
public String toString(final Object o);
}
