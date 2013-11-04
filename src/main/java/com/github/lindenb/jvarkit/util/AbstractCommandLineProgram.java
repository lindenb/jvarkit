package com.github.lindenb.jvarkit.util;

import java.io.InputStream;
import java.net.InetAddress;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.jar.Manifest;
import java.util.logging.Logger;

public abstract class AbstractCommandLineProgram
	{
	private static final Logger LOG=Logger.getLogger("jvarkit");
	private String commandLine="";
	private String version=null;
	private String compileDate;
	
	
	protected AbstractCommandLineProgram()
		{
		}
	
	protected void printStandardPreamble()
		{
		System.out.println(getClass().getSimpleName());
		System.out.println("Compilation "+getCompileDate());
		System.out.println("Version "+getVersion());
		System.out.println(getProgramDescription());
		
		}
	public void printUsage()
		{
		printStandardPreamble();
		System.out.println( "See https://github.com/lindenb/jvarkit");
		}
	
	public String getProgramDescription()
		{
		return getClass().getName();
		}
	
    
	public Logger getLogger()
		{
		return LOG;
		}
	
	public abstract int doWork(String args[]);
	
	private void loadManifest()
		{
		try
			{
			
				Enumeration<URL> resources = getClass().getClassLoader()
						  .getResources("META-INF/MANIFEST.MF");//not '/META-INF'
				while (resources.hasMoreElements())
					{
					URL url=resources.nextElement();
					InputStream in=url.openStream();
					if(in==null)
						{
						continue;
						}
					
					Manifest m=new Manifest(in);
					in.close();
					in=null;
					java.util.jar.Attributes attrs=m.getMainAttributes();
					if(attrs==null)
						{
						continue;
						}
					String s =attrs.getValue("Git-Hash");
					if(s!=null && !s.isEmpty() && !s.contains("$")) //ant failed
						{
						this.version=s;
						}
					
					s =attrs.getValue("Compile-Date");
					if(s!=null && !s.isEmpty()) //ant failed
						{
						this.compileDate=s;
						}
					}
			}	
		catch(Exception err)
			{
			
			}
		
		}
	
	public String getCompileDate()
		{
		if(this.compileDate==null)
			{
			this.compileDate="undefined";
			loadManifest();
			}
		return compileDate;
		}
	
	public String getVersion()
		{
		if(this.version==null)
			{
			this.version="1.0";
			loadManifest();
			}
		return version;
		}
 
	protected int instanceMain(String args[])
		{
		StringBuilder b=new StringBuilder();
		for(int i=0;i< args.length;++i)
			{
			if(i!=0) b.append(' ');
			b.append(args[i]);
			}
		this.commandLine=b.toString();
		b=null;
		Date startDate=new Date();
		getLogger().info("Starting JOB at "+startDate+" "+getClass().getName()+
				" version="+getVersion()+" "+
				" built="+getCompileDate());
		getLogger().info(this.commandLine);
		String hostname="";
		try
			{
			hostname= InetAddress.getLocalHost().getHostName();
			}
		catch(Exception err)
			{
			hostname="host";
			}
		getLogger().info("Executing as " +
                System.getProperty("user.name") + "@" + hostname +
                " on " + System.getProperty("os.name") + " " + System.getProperty("os.version") +
                " " + System.getProperty("os.arch") + "; " + System.getProperty("java.vm.name") +
                " " + System.getProperty("java.runtime.version") 
                );
		int ret=doWork(args);
	
		final Date endDate = new Date();
		final double elapsedMinutes = (endDate.getTime() - startDate.getTime()) / (1000d * 60d);
        final String elapsedString  = new DecimalFormat("#,##0.00").format(elapsedMinutes);
    	getLogger().info("End JOB status="+ret+" [" + endDate + "] " + getClass().getName() + " done. Elapsed time: " + elapsedString + " minutes.");
		return ret;
		}
	
    public void instanceMainWithExit(final String[] argv) {
        System.exit(instanceMain(argv));
    }
}
