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
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC
## Example

```bash
$   java -jar dist/vcfjmx.jar -p MyWorkflow1 input.vcf > /dev/null
```

while the stream is running, open a new jconsole https://docs.oracle.com/javase/7/docs/technotes/guides/management/jconsole.html . here you can get the number of records, the elapsed time. Two operation are available:

* doBreak: interrupt current streaming , exit with success (0)
* doAbort: interrupt current streaming , exit with failure (-1)

```
$ java -jar dist/vcfjmx.jar -p 1000G  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz | \
    java -jar dist/vcfjmx.jar  -p 1000G-2  |\
    java -jar dist/vcfjmx.jar  -p 1000G-3  |\
    java -jar dist/vcfjmx.jar  -p 1000G-4 > /dev/null

[INFO/VcfJmx] 2015-07-10 14:11:08 "Starting JOB at Fri Jul 10 14:11:08 CEST 2015 com.github.lindenb.jvarkit.tools.jmx.VcfJmx version=4f797a9fbf2c3ceac9cec3c431c719ad794953c2  built=2015-07-10:13-07-05"
[INFO/VcfJmx] [INFO/VcfJmx] 2015-07-10 14:11:08 "Command Line args : -p 1000G-2"
2015-07-10 14:11:08 "Starting JOB at Fri Jul 10 14:11:08 CEST 2015 com.github.lindenb.jvarkit.tools.jmx.VcfJmx version=4f797a9fbf2c3ceac9cec3c431c719ad794953c2  built=2015-07-10:13-07-05"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Command Line args : -p 1000G-3"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Starting JOB at Fri Jul 10 14:11:08 CEST 2015 com.github.lindenb.jvarkit.tools.jmx.VcfJmx version=4f797a9fbf2c3ceac9cec3c431c719ad794953c2  built=2015-07-10:13-07-05"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Command Line args : -p 1000G ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19"
[INFO/VcfJmx] 2015-07-10 14:11:08 "Executing as lindenb@kaamelot-master01 on Linux 2.6.32-431.17.1.el6.x86_64 amd64; Java HotSpot(TM) 64-Bit Server VM 1.7.0_60-b19"
[INFO/VcfJmx] 2015-07-10 14:11:08 "reading from stdin"
[INFO/VcfJmx] 2015-07-10 14:11:08 "reading from stdin"
[INFO/VcfJmx] 2015-07-10 14:11:08 "reading from ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz"

[SEVERE/VcfJmx] 2015-07-10 14:11:46 "#### Process "1000G-3" received message BREAK"

[INFO/VcfJmx] 2015-07-10 14:11:46 "Number of Variants:1130774"
[INFO/VcfJmx] 2015-07-10 14:11:46 "End JOB status=0 [Fri Jul 10 14:11:46 CEST 2015] com.github.lindenb.jvarkit.tools.jmx.VcfJmx done. Elapsed time: 0.64 minutes."
[INFO/VcfJmx] 2015-07-10 14:11:46 "Number of Variants:1120774"
[INFO/VcfJmx] 2015-07-10 14:11:46 "Number of Variants:1110975"
[INFO/VcfJmx] 2015-07-10 14:11:46 "End JOB status=0 [Fri Jul 10 14:11:46 CEST 2015] com.github.lindenb.jvarkit.tools.jmx.VcfJmx done. Elapsed time: 0.64 minutes."
[INFO/VcfJmx] 2015-07-10 14:11:46 "End JOB status=0 [Fri Jul 10 14:11:46 CEST 2015] com.github.lindenb.jvarkit.tools.jmx.VcfJmx done. Elapsed time: 0.64 minutes."
```
END_DOC

 */
@Program(name="vcfjmx",
description="Monitor/interrupt/break a VCF stream with java JMX http://www.oracle.com/technetwork/articles/java/javamanagement-140525.html",
keywords={"java","jmx","vcf"})
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
	protected int doVcfToVcf(String inputName, VCFIterator in,
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
		new VcfJmx().instanceMainWithExit(args);
		}

}
