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

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class TeeInputStream extends InputStream
		{
		private  InputStream in;
	    private  OutputStream out;
	    /** shall we close out after 'in' ? */
	    private boolean closeOutOnCloseIn=false;

	  
	    public TeeInputStream(final InputStream in,final OutputStream out, boolean closeOutOnCloseIn)
	    	{
	        this.in = in;
	        this.out=out;
	        this.closeOutOnCloseIn = closeOutOnCloseIn;
	    	}

	    @Override
	    public void close() throws IOException
	    	{
	        try {
	            in.close();
	        	} finally {
	            if (closeOutOnCloseIn && out!=null)
	            	{
	            	out.flush();
	            	out.close();
	            	}
	        	}
	    	}

	    public int read() throws IOException
	    	{
	        int ch = this.in.read();
	        if (ch != -1 && out!=null) {
	        	 out.write(ch);
	        	}
	        return ch;
	    	}

	
	    public int read(byte[] bts, int st, int end) throws IOException {
	        int n = super.read(bts, st, end);
	        if (n != -1 && out!=null)
	        	{
	            out.write(bts, st, n);
	        	}	
	        return n;
	    }

	   
	    public final int read(byte[] bts) throws IOException
	    	{
	        return read(bts,0,bts.length);
	    	}

		}
