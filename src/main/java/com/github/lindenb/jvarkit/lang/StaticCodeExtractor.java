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
package com.github.lindenb.jvarkit.lang;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Optional;

import com.github.lindenb.jvarkit.util.log.Logger;


public class StaticCodeExtractor {
	private static final Logger LOG = Logger.build(StaticCodeExtractor.class).make();

	private Class<?>  clazz;
	private StaticCodeExtractor() {
		}
	public static StaticCodeExtractor forClass(final Class<?> clazz) {
		if(clazz==null) throw new IllegalArgumentException("class is null");
		final StaticCodeExtractor ex = new StaticCodeExtractor();
		ex.clazz = clazz;
		return ex;
		}
	public Optional<String> extract(final String opcode) {
		if(StringUtils.isBlank(opcode)) throw new IllegalArgumentException("opcode is null or empty");
		final String begin = "BEGIN_"+opcode;
		final String end = "END_"+opcode;
		String className=  this.clazz.getName();
		final int dollar=className.indexOf('$');
		if(dollar!=-1) className=className.substring(0, dollar);
		className=className.replace('.', '/')+".java";

		InputStream in= clazz.getResourceAsStream("/"+className);
		if(in==null) return Optional.empty();
		boolean gotit = false;
		try {
			StringBuilder sb = null;
			try(BufferedReader br= new BufferedReader(new InputStreamReader(in))) {
				String line;
				
				while((line=br.readLine())!=null) {
					if(line.equals(begin)) {
						sb = new StringBuilder();
						while((line=br.readLine())!=null) {
							if(line.equals(end)) {
								gotit=true;	
								break;
								}
							sb.append(line).append("\n");
							}
						break;
						}
					}
				}
			if(!gotit) return Optional.empty();
			return Optional.of(sb.toString());
			} 
		catch(final IOException err) {
			LOG.error(err);
			return Optional.empty();	
			}
		finally {
			try {in.close();} catch(Throwable err) {}
			}
		}
	}
