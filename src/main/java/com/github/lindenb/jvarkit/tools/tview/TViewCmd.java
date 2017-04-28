/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.tview;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;


public class TViewCmd extends Launcher
	{
	private static final Logger LOG = Logger.build(TViewCmd.class).make();
	@Parameter(names={"-o","--out"},description="Output file . Default: stdout")
	private File outputFile = null;
	@Parameter(names={"-r","--region"},description="Interval list")
	private Set<String> intervalStr = new HashSet<>();
	@Parameter(names={"-width","--width"},description="default screen width")
	private int defaultScreenWidth=-1;
	
	@ParametersDelegate
	TView tview = new TView();

	private  Interval parseInterval(String regionStr) {
		final int colon=regionStr.indexOf(':');
		;
		int chromStart;
		int chromEnd;

		if(colon==-1) throw new IllegalArgumentException("colon missing in "+regionStr);
		final String chrom=regionStr.substring(0,colon);
			
		
		final int hyphen=regionStr.indexOf('-', colon+1);
		if(hyphen==-1)
			{
			int cols = defaultScreenWidth;
			
			chromStart=Integer.parseInt(regionStr.substring(colon+1));
			String COLUMNS=System.getenv("COLUMNS");
			if(COLUMNS!=null && COLUMNS.matches("[0-9]+")) {
				cols = Integer.parseInt(COLUMNS);
				LOG.info(cols);
				}
			if(cols<=1) cols=80;
			chromEnd = chromStart + cols;
			}
		else
			{
		
			chromStart=Integer.parseInt(regionStr.substring(colon+1,hyphen));
			chromEnd=Integer.parseInt(regionStr.substring(hyphen+1));
			}
		if(chromStart<0 || chromEnd<chromStart)
			{
			throw new IllegalArgumentException("bad position in "+regionStr);
			}
		return new Interval(chrom, chromStart, chromEnd,false,regionStr);
		}
	
	@Override
	public int doWork(List<String> args) {
		PrintStream out = null;
		try
			{	
			args = new ArrayList<>(IOUtils.unrollFiles(args));
			if(args.isEmpty())
				{
				LOG.error("no BAM file defined");
				return -1;
				}
			tview.setBamFiles(args.stream().map(S->SamInputResource.of(S)).
					collect(Collectors.toList()));
			if(tview.initialize()!=0)
				{
				LOG.error("cannot initialize tview");
				return -1;
				}
			out = super.openFileOrStdoutAsPrintStream(outputFile);
			for(final Interval interval : this.intervalStr.stream().map(S->parseInterval(S) ).collect(Collectors.toList())) {
				tview.setInterval(interval);
				tview.paint(out);
			}
			
			out.flush();
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(tview);
			}
		
		}
	
	public static void main(String[] args)
		{
		new TViewCmd().instanceMainWithExit(args);
		}
	}
