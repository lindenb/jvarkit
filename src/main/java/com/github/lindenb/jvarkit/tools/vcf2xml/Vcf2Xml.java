/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.io.File;
import java.io.IOException;
import java.util.List;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.XMLVcfWriterFactory;
/**

 
 */
@Program(name="vcf2xml",description="Convert VCF to XML",keywords={"vcf","xml"})
public class Vcf2Xml extends Launcher
	{

	@Parameter(names={"-o","--output"},description="Output file. Optional . Default: stdout")
	private File outputFile = null;

	public Vcf2Xml()
		{
		
		}
	@Override
	protected int doVcfToVcf(
			final String inputName,
			final VcfIterator iterin, 
			final VariantContextWriter out
			) {
		final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(iterin.getHeader());
		while(iterin.hasNext())
			{
			out.add(progress.watch(iterin.next()));
			}
		progress.finish();
		return 0;
		}
	
	/** open VariantContextWriter */
	@Override
	protected VariantContextWriter openVariantContextWriter(final File outorNull) throws IOException {
		final XMLVcfWriterFactory factory=XMLVcfWriterFactory.newInstance();
		if(outorNull!=null)
			{
			factory.setOutputFile(outorNull);
			}
		return factory.createVariantContextWriter();
		}
	
	@Override
	public int doWork(List<String> args) {
		return doVcfToVcf(args, outputFile);
		}
	
	public static void main(final String[] args)
		{
		new Vcf2Xml().instanceMain(args);
		}
}
