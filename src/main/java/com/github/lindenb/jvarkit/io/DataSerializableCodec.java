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
package com.github.lindenb.jvarkit.io;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.util.function.Supplier;

import com.github.lindenb.jvarkit.util.picard.AbstractDataCodec;


/** interface for object I/O to/from DataInput/OutputStream */
public class DataSerializableCodec<T extends DataSerializable> 
	extends AbstractDataCodec<T>{
	private Supplier<T> factory;
	public DataSerializableCodec(final Supplier<T> factory) {
		this.factory = factory;
		}
	public DataSerializableCodec(final Class<T> clazz) {
		try{
			final Constructor<T> ctor = clazz.getConstructor();
			this.factory=()->{
				try {
					return ctor.newInstance();
					}
				catch(final Exception err)
					{
					throw new RuntimeException(err);
					}
				};
			}
		catch(final Exception err)
			{
			throw new RuntimeException(err);
			}
		}
	public Supplier<T> getFactory() {
		return factory;
		}
	@Override
	public void encode(final DataOutputStream dos, final T object) throws IOException {
		object.writeDataOutputStream(dos);
		}
	
	@Override
	public T decode(final DataInputStream dis) throws IOException {
		try
			{
			final T o = getFactory().get();
			o.readDataInputStream(dis);
			return o;
			}
		catch(final EOFException err) {
			return null;
			}
		}
	
	@Override
	public DataSerializableCodec<T> clone() {
		return new DataSerializableCodec<>(getFactory());
		}
	}
