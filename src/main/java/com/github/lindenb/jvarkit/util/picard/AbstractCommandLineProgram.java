package com.github.lindenb.jvarkit.util.picard;

import java.io.IOException;
import java.io.InputStream;
import java.util.jar.Manifest;


import net.sf.picard.cmdline.CommandLineProgram;

public abstract class AbstractCommandLineProgram
	extends CommandLineProgram
	{
	@Override
	public String getVersion()
		{
		return getGitHash();
		}
	@Override
	public String getProgramVersion()
		{
		return getGitHash();
		}

	public String  getGitHash()
		{
		String hash="";
		InputStream in=getClass().getResourceAsStream("META-INF/MANIFEST.MF");
		if(in!=null)
			{
			try
				{
				Manifest m=new Manifest(in);
				java.util.jar.Attributes attrs =m.getMainAttributes();
				hash=attrs.getValue("Git-Hash");
				if(hash==null) hash="";
				if(hash.contains("$")) //ant failed
					{
					hash="";
					}
				}
			catch(Exception err)
				{
				hash="";
				}
			finally
				{
				try { in.close();}
				catch(IOException err) {}
				}
			}
		return hash;
		}
	
	@Override
	protected boolean parseArgs(final String[] argv)
		{
		return true;
		}
	
	
	}
