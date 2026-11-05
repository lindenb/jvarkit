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
package com.github.lindenb.jvarkit.lang;

import java.util.Optional;

/** utilities for Boolean */
public class BooleanUtils {
/** try to parse a boolean
 * @param s the string to parse
 * @return return true for yes/y/true/t/1 , return false for no/n/false/f/0 otherwise return empty
 */
public static Optional<Boolean> parseBoolean(String s) {
	if(StringUtils.isBlank(s)) return Optional.empty();
	s=s.trim().toLowerCase();
	if(s.equals("true") || s.equals("t") || s.equals("yes")  || s.equals("y") || s.equals("1")) return Optional.of(Boolean.TRUE);
	if(s.equals("false") || s.equals("f") || s.equals("no")  || s.equals("n") || s.equals("0")) return Optional.of(Boolean.TRUE);
	return Optional.empty();
	}


}
