package com.github.lindenb.jvarkit.tools.cgi;

import java.io.IOException;
import java.io.InputStream;

import htsjdk.samtools.util.CloserUtil;

public abstract class AbstractCGICallApp extends AbstractCGI
	{
	protected static class ConsummeInputStreamThread extends Thread	
	   {
	   private InputStream is;
	   public ConsummeInputStreamThread(InputStream is)
	    	{
	        this.is = is;
	    	}
	  
	    @Override
	    public void run()
		    {
	        try
		        {
	        	while (is.read()!=-1);
		        }
	        catch (IOException ioe)
              	{
        		
              	}
	        finally
	        	{
	        	CloserUtil.close(this.is);
	        	}
		    }
	    
		}
	
	protected AbstractCGICallApp()
		{
		}

	}
