/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tests;

import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;


public class TestRunnerImpl implements TestRunner {
	private final Set<Class<?>> classToTest = new HashSet<>();
	private int level=0;
	private int nClasses=0;
	private int nMethods=0;
	private int nSuccess=0;
	private int nFailures=0;
	
	private void begin(final String title) {
		
		++level;
	}
	
	private void end() {
		--level;
	}
	
	
	public void register(final Class<?> c)  {
		if(c==null || this.classToTest.contains(c)) return;
		this.classToTest.add(c);
		final AlsoTest also = c.getDeclaredAnnotation(AlsoTest.class);
		if(also!=null) {
			for(Class<?> c2:also.value()) {
				register(c2);
				}
			}
		Class<?> parent = c.getSuperclass();
		if(parent!=null && !parent.equals(Object.class)) {
			register(parent);
			}
		}
	
	private void run(Method m) {
		this.nMethods++;
		begin(m.getName());
		try {
			m.invoke(null, this);
			this.nSuccess++;
			}
		catch(final Throwable err) {
			
			this.nFailures++;
			}
		end();
		}
	
	private void run(Class<?> C) {
		this.nClasses++;
		begin(C.getName());
		for(Method m:C.getDeclaredMethods()) {
			final Class<?> paramTypes[] = m.getParameterTypes();
			if(!m.getName().startsWith("test")) continue;
			if(paramTypes.length!=1) continue;
			if(paramTypes[0].equals(TestRunner.class)) continue;
			if(!Modifier.isStatic(m.getModifiers())) continue;
			m.setAccessible(true);
			run(m);
		}
		end();
	}
	
	@Override
	public void run(final Path path) {
		this.nClasses = 0;
		this.nMethods = 0;
		begin(String.valueOf(this.classToTest.size())+" classes(s)");
		for(final Class<?> c: this.classToTest) {
			run(c);
		}
		end();
	}
}
