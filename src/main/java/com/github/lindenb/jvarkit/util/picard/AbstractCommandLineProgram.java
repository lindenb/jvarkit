package com.github.lindenb.jvarkit.util.picard;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.jar.Manifest;


import htsjdk.samtools.cmdline.CommandLineProgram;

public abstract class AbstractCommandLineProgram
	extends CommandLineProgram
	{
	
	@Override
	public String getVersion()
		{
		String s=getGitHash();
		return (s.isEmpty()?"null":s);
		}
	@Override
	public String getProgramVersion()
		{
		String s=getGitHash();
		return (s.isEmpty()?"null":s);
		}

	public String  getGitHash()
		{
		String hash="";
		InputStream in=getClass().getResourceAsStream("/META-INF/MANIFEST.MF");
		if(in==null) return hash;
			
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
			
		return hash;
		}
	
	protected void  testRemoteGit()
		{
		String hash=getGitHash();
		if(hash.isEmpty()) return;
		hash="\""+hash+"\"";
		InputStream in=null;
		try
			{
			URL url=new URL("https://api.github.com/repos/lindenb/jvarkit/git/refs/heads/master");
			URLConnection con=url.openConnection();
			con.setReadTimeout(3*1000);
			con.connect();
			in=con.getInputStream();
			StringBuilder b=new StringBuilder();
			int c;
			while((c=in.read())!=-1) b.append((char)c);
			if(b.toString().indexOf(hash)==-1)
				{
				System.err.println("This program (hash:"+hash+") is not synchronized with the git repository : "+b);
				return;
				}
			}
		catch(Exception err)
			{
			
			}
		finally
			{
			try {if(in!=null) in.close();}
			catch(IOException err) {}
			}
		
		}
	
	
	}
