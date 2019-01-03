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
package com.github.lindenb.jvarkit.util.ncbi;


/**
 * Constants for the NCBI.
 */
public class NcbiConstants {
	public static final String EUTILS_BASE_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";

	private static String service(final String name) {
		return EUTILS_BASE_URL+name+".fcgi";
		}
	/** return service base url */
	public static String efetch() {
		return service("efetch");
		}
	/** return service base url */
	public static String esearch() {
		return service("esearch");
		}
	/** return service base url */
	public static String elink() {
		return service("elink");
		}
	/** return service base url */
	public static String esummary() {
		return service("esummary");
		}
	/** 
	 * " Increasing retmax allows more of the retrieved UIDs to be included in the XML output, up to a maximum of 100,000 records" 
	 * https://www.ncbi.nlm.nih.gov/books/NBK25499/
	 * */
	public static final int RETMAX_MAX = 100_000;
	
}
