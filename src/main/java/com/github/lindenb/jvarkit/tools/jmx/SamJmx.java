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
package com.github.lindenb.jvarkit.tools.jmx;

import java.io.File;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.util.List;

import javax.management.MBeanServer;
import javax.management.ObjectName;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloserUtil;

/*
BEGIN_DOC

## Example

```bash
$   java -jar dist/samjmx.jar  -T bam -p MyWorkflow1 input.bam > /dev/null
```

while the stream is running, open a new jconsole https://docs.oracle.com/javase/7/docs/technotes/guides/management/jconsole.html . here you can get the number of records, t
he elapsed time. Two operation are available:

* doBreak: interrupt current streaming , exit with success (0)
* doAbort: interrupt current streaming , exit with failure (-1)



END_DOC
*/
@Program(name="samjmx",
description="Monitor/interrupt/break a BAM/SAM stream with java JMX http://www.oracle.com/technetwork/articles/java/javamanagement-140525.html",
keywords={"sam","bam","jmx","monitoring"})
public class SamJmx extends Launcher
	{
	private static final Logger LOG=Logger.build(SamJmx.class).make();
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File output=null;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs = new WritingBamArgs();
	@Parameter(names={"-p"},required=false,description="Stream identifier")
	private String projectName=null;
	
	public SamJmx() {
	}
	
	public void setProjectName(final String projectName) {
		this.projectName = projectName;
		}
	
	private int doWork(final SamReader in) throws IOException
		{
		String name=this.projectName;
		if(name==null || name.trim().isEmpty())
			{
			name= "undefined";
			}
		final MBeanServer mbeanServer = ManagementFactory.getPlatformMBeanServer();
		final LocatableStreamInfo dynamicMBean=new LocatableStreamInfo( name );
		LocatableStreamInfo.Status status = LocatableStreamInfo.Status.IDLE;
		ObjectName objectMBean = null;
		SAMFileWriter out=null;
		SAMRecordIterator iter=null;
		try
			{
		    out = this.writingBamArgs.openSAMFileWriter(output,in.getFileHeader(), true);
		    
		    objectMBean = new ObjectName(
		    		dynamicMBean.getClass().getPackage()
		            .getName() + ":type=" + dynamicMBean.getClass().getSimpleName()
		            
		    		);
		    mbeanServer.registerMBean(dynamicMBean,objectMBean);
			iter = in.iterator();
			while(iter.hasNext())
				{
				final SAMRecord rec = iter.next();
				out.addAlignment(rec);
				status = dynamicMBean.watch(rec);
				if( status == LocatableStreamInfo.Status.ABORT ||
				    status == LocatableStreamInfo.Status.BREAK)
					{
					LOG.error("#### Process \""+name+"\" received message "+status.name());
					break;
					}
				}
			if(status == LocatableStreamInfo.Status.ABORT)
				{
				LOG.error("#### Process \""+name+"\" : Exit failure");
				mbeanServer.unregisterMBean(objectMBean);
				if(this.output!=null)
					{
					CloserUtil.close(out);
					out=null;
					this.output.delete();
					}
				System.exit(-1);
				}
			return 0;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(iter);
			CloserUtil.close(out);
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
		SamReader in=null;
		try
			{
			in = super.openSamReader(super.oneFileOrNull(args));
			return doWork(in); 
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			}
		}

	public static void main(final  String[] args) throws IOException
		{
		new SamJmx().instanceMainWithExit(args);
		}
	}
