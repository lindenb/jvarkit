/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.util.Collection;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.util.command.Command;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.XMLVcfWriterFactory;
/**
 * @author lindenb
 *
 */
public class Vcf2Xml extends AbstractVcf2Xml
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(Vcf2Xml.class);

	@Override
	public Command createCommand() {
		return new MyCommand();
		}

	static private class MyCommand extends AbstractVcf2Xml.AbstractVcf2XmlCommand
		{    
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
		{
		VcfIterator in=null;
		VariantContextWriter w = null;
		try {
			in = openVcfIterator(inputName);
			XMLVcfWriterFactory factory=XMLVcfWriterFactory.newInstance();
			if(getOutputFile()!=null)
				{
				w = factory.createVariantContextWriter(getOutputFile());
				}
			else
				{
				w = factory.createVariantContextWriter(stdout(),"UTF-8");
				}
			w.writeHeader(in.getHeader());
			SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(in.getHeader());
			while(in.hasNext())
				{
				w.add(progress.watch(in.next()));
				}
			progress.finish();
			LOG.info("done");
			return RETURN_OK;
		} catch (Exception e) {
			e.printStackTrace();
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(w);
			}
		}
	}
	
	public static void main(String[] args)
		{
		new Vcf2Xml().instanceMain(args);
		}
	}
