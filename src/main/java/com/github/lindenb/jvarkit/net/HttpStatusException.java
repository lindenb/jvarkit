/*
The MIT License (MIT)

Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.net;

import org.apache.http.HttpException;
import org.apache.http.HttpStatus;
import org.apache.http.StatusLine;

/**
 * 
 * HttpStatusException
 *
 */
@SuppressWarnings("serial")
public class HttpStatusException extends HttpException {
	private static final long serialVersionUID = 0L;
	private final StatusLine statusLine;
	public HttpStatusException(final StatusLine statusLine) {
		super("Error "+statusLine.getStatusCode()+":"+statusLine.getReasonPhrase());
		this.statusLine = statusLine;
		}
	public StatusLine getStatusLine() {
		return statusLine;
		}
	public int getStatusCode() {
		return getStatusLine().getStatusCode();
		}
	
	public StatusLine assertHttpStatusOK(StatusLine sl) throws HttpStatusException {
		if(sl.getStatusCode()!=HttpStatus.SC_OK) throw new HttpStatusException(sl);
		return sl;
		}
	}
