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
package com.github.lindenb.jvarkit.util.vcf;


/**
 * Author: Pierre Lindenbaum PhD. @yokofakun
 * Motivation http://www.biostars.org/p/66319/ 
 */

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.util.AbstractCommandLineProgram;


/**
 * extends com.github.lindenb.jvarkit.util.AbstractCommandLineProgram
 * while AbstractVCFFilter extends com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram
 * @author lindenb
 *
 */
@Deprecated
public abstract class AbstractVCFFilter2
	//extends AbstractCommandLineProgram
	{
	
	private AbstractVCFFilter2()
		{
		
		}
	
	public String getProgramDescription() {
		return "Another VCF filter.";
		}
	
	protected abstract void doWork(
			VcfIterator in,
			VariantContextWriter out
			) throws IOException;
	
	protected  VcfIterator createVcfIterator(String IN) throws IOException
		{
		return VCFUtils.createVcfIterator(IN);
		}
	
	protected VariantContextWriter createVariantContextWriter(File OUT) throws IOException
		{
		return VCFUtils.createVariantContextWriter(OUT);
		}
	
	
	/** returns wether the output (stdout) has raised an error */
	protected boolean outCheckError()
		{
		return System.out.checkError();
		}
	
	protected int doWork(String IN,File OUT)
		{
		VcfIterator r=null;
		VariantContextWriter w=null;
		try
			{
			r=this.createVcfIterator(IN);
			w=this.createVariantContextWriter(OUT);
			doWork(r,w);
			}
		catch (Exception e)
			{
			return -1;
			}
		finally
			{
			CloserUtil.close(w);
			CloserUtil.close(r);
			}	
		return 0;
		}
	
	protected int doWork(int optind,String args[])
		{
		return -1;
		}
	
	}
