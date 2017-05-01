/*
The MIT License (MIT)

Copyright (c) 2015 Pierre Lindenbaum

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


History:
* 2015 creation

*/
package com.github.lindenb.jvarkit.tools.jmx;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.List;

import javax.management.MBeanServer;
import javax.management.ObjectName;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

@Program(name="vcfjmx",description="Monitor/interrupt/break a VCF stream with java JMX http://www.oracle.com/technetwork/articles/java/javamanagement-140525.html")
public class VcfJmx extends Launcher
	{
	private static final Logger LOG=Logger.build(VcfJmx.class).make();
	@Parameter(names={"-o","--out"},required=false,description="Output vcf , ot stdin")
	private File outputFile=null;
	@Parameter(names={"-p"},required=false,description="Stream identifier")
	private String projectName=null;
	
	public VcfJmx() {
	}
	
	public void setProjectName(String projectName) {
		this.projectName = projectName;
		}
	
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in,
			VariantContextWriter out)
		{
		String name=this.projectName;
		if(name==null || name.trim().isEmpty())
			{
			name= "undefined";
			}
	    MBeanServer mbeanServer = ManagementFactory.getPlatformMBeanServer();
		LocatableStreamInfo dynamicMBean=new LocatableStreamInfo( name );
		LocatableStreamInfo.Status status = LocatableStreamInfo.Status.IDLE;
		ObjectName objectMBean = null;
		try
			{
		    
		  
		    objectMBean = new ObjectName(
		    		dynamicMBean.getClass().getPackage()
		            .getName() + ":type=" + dynamicMBean.getClass().getSimpleName()
		            
		    		);
		    mbeanServer.registerMBean(dynamicMBean,objectMBean);
			out.writeHeader(in.getHeader());
			while(in.hasNext())
				{
				VariantContext ctx = in.next();
				out.add(ctx);
				status = dynamicMBean.watch(ctx);
				if( status == LocatableStreamInfo.Status.ABORT ||
				    status == LocatableStreamInfo.Status.BREAK)
					{
					LOG.error("#### Process \""+name+"\" received message "+status.name());
					break;
					}
				if(System.out.checkError()) break;
				}
			if(status == LocatableStreamInfo.Status.ABORT)
				{
				LOG.error("#### Process \""+name+"\" : Exit failure");
				mbeanServer.unregisterMBean(objectMBean);
				System.exit(-1);
				}
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			if(objectMBean!=null)
				{
				try { mbeanServer.unregisterMBean(objectMBean);}
				catch(Exception err2) {}
				}
			}
		}

	@Override
	public int doWork(final List<String> args)
		{
		return doVcfToVcf(args, outputFile);
		}

	
	public static void main(String[] args) throws IOException
		{
		new VcfJmx().instanceMain(args);
		}

}
