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

	  
	    public TeeInputStream( InputStream in, OutputStream out, boolean closeOutOnCloseIn)
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
