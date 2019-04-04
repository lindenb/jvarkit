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
package com.github.lindenb.jvarkit.tools.misc;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.springframework.context.ApplicationContext;
import org.springframework.context.support.FileSystemXmlApplicationContext;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
/**
BEGIN_DOC





END_DOC

*/
@Program(name="vcfspringfilter",
	description="Uses the java spring Framework ( https://docs.spring.io/spring/docs/current/spring-framework-reference/htmlsingle/ ) to build complex vcf filters",
	keywords={"vcf","java","spring","framework"},
	generate_doc=false
	)
public class VcfSpringFilter extends Launcher {
	private static final Logger LOG=Logger.build(VcfSpringFilter.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile=null;
	@Parameter(names={"-c","--config"},description="Spring configuration file.")
	private List<String> springCongigFiles = new ArrayList<>();
	@Parameter(names={"-m","--main"},description="Main bean name")
	private String mainBeanName="main";

	
	interface VcfChain
		{
		int doFilter(VCFIterator in,VariantContextWriter out);
		}
	
	private ApplicationContext springApplicationContext = null;
	
	public VcfSpringFilter()
		{
		
		}
	@Override
	public int doWork(final List<String> args) {
		VCFIterator vcfIn  = null;
		VariantContextWriter vcfOut=null;
		try {
			if(this.springCongigFiles.isEmpty())
				{
				LOG.error("no spring config file was provided");
				return -1;
				}
			this.springApplicationContext = 
					new FileSystemXmlApplicationContext(
							springCongigFiles.toArray(new String[springCongigFiles.size()])
					);

			
			if(!this.springApplicationContext.containsBean(this.mainBeanName))
				{
				LOG.error("cannot get bean "+mainBeanName+" in "+
						this.springCongigFiles
						); 
				return -1;
				}
			final Object o = this.springApplicationContext.getBean(mainBeanName);
			if( o == null || !(o instanceof VariantContextWriter))
				{
				LOG.error("bean "+mainBeanName+" is not a  VcfChain but "+
						(o==null?"null":o.getClass().getName())
						); 
				return -1;
				}
			final VariantContextWriter writer = VariantContextWriter.class.cast(o);
			
			final String inputFile= oneFileOrNull(args);
			vcfIn = super.openVCFIterator(inputFile);
			
			
			return ret;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(vcfIn);
			CloserUtil.close(vcfOut);
			this.springApplicationContext=null;
			}
		}
	
	public static void main(final String[] args) {
		new VcfSpringFilter().instanceMainWithExit(args);

	}

}
