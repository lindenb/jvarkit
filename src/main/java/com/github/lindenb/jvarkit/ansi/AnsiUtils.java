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
package com.github.lindenb.jvarkit.ansi;

import com.github.lindenb.jvarkit.lang.StringUtils;

public class AnsiUtils {
	public static final String ANSI_ESCAPE = "\u001B[";
	public static final String ANSI_RESET = ANSI_ESCAPE+"0m";
	public static char HISTOGRAM_CHARS[] = new char[]{
			' ',
			'\u2581', '\u2582', '\u2583', '\u2584', '\u2585', 
			'\u2586', '\u2587', '\u2588'
			};
	
		public static enum AnsiColor {
			BLACK (30),
			RED (31),
			GREEN (32),
			YELLOW (33),
			BLUE (34),
			MAGENTA (35),
			CYAN (36),
			WHITE (37)
			;
		
		AnsiColor(final int opcode) {
			this.opcode=opcode;
			}
		final int opcode;
		
		public String begin() {
			return ANSI_ESCAPE+this.opcode+"m";
			}
		public String end() {
			return ANSI_RESET;
			}
		public String colorize(final String str) {
			if(StringUtils.isBlank(str)) return str;
			return begin() + str + end();
			}
		}
	
	
	public static char getHistogram(double fraction) {
		if(fraction<=0) return HISTOGRAM_CHARS[0];
		switch((int)(fraction*10.0)) {
			case 0: return HISTOGRAM_CHARS[0];
			case 1:  case 2: return HISTOGRAM_CHARS[1];
			case 3: return HISTOGRAM_CHARS[2];
			case 4: return HISTOGRAM_CHARS[3];
			case 5: return HISTOGRAM_CHARS[4];
			case 6: return HISTOGRAM_CHARS[5];
			case 7: return HISTOGRAM_CHARS[6];
			case 8: return HISTOGRAM_CHARS[7];
			default: return HISTOGRAM_CHARS[8];
			}
		}
	}
