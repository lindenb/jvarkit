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
package com.github.lindenb.jvarkit.jexl;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.function.Function;

import org.apache.commons.jexl2.JexlContext;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

/**
 * A JEXL Expression handler returning a string
 */
public class JexlToString extends BaseJexlHandler
	implements Function<JexlContext,String>
	{
	/** compile JexlExpressionToString as {@link org.apache.commons.jexl2.Expression} **/
	public JexlToString(final String expr) {
		this(expr,true);
		}
	
	public JexlToString(final String expr,boolean is_expression) {
		super(expr,is_expression);
		}
	/** output is never null, null is converted to empty string, java.lang.Object is converted with toString() */
	@Override
	public String apply(final JexlContext t) {
		final Object o = eval(t);
		if(o==null) return "";
		return o.toString();
		}
	
	/** load script stored in path and compile it */
	public static JexlToString compileScript(final Path path) {
		IOUtil.assertFileIsReadable(path);
		try {
			return new JexlToString(Files.readString(path),false);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		}
	}
