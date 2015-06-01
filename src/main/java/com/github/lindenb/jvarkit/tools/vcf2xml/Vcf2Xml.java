/**
 * 
 */
package com.github.lindenb.jvarkit.tools.vcf2xml;

import java.io.IOException;
import java.io.PrintStream;

import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter3;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.XMLVcfWriterFactory;
/**
 * @author lindenb
 *
 */
public class Vcf2Xml extends AbstractVCFFilter3
	{
	public Vcf2Xml()
		{
		
		}
	@Override
	protected void doWork(String inputSource, VcfIterator vcfIn,
		VariantContextWriter out) throws IOException
		{
		SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(vcfIn.getHeader());
		while(vcfIn.hasNext())
			{
			out.add(progress.watch(vcfIn.next()));
			}
		progress.finish();
		}
	
	/** open VariantContextWriter */
	@Override
	protected VariantContextWriter createVariantContextWriter()
		throws IOException
		{
		XMLVcfWriterFactory factory=XMLVcfWriterFactory.newInstance();
		if(getOutputFile()!=null)
			{
			factory.setOutputFile(getOutputFile());
			}
		return factory.createVariantContextWriter();
		}

	
	
	@Override
	public String getProgramDescription() {
		return "Convert VCF to XML";
		}
	
	@Override
    protected String getOnlineDocUrl() {
    	return DEFAULT_WIKI_PREFIX+"Vcf2Xml";
    }
	
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -o (out)  output file. default stdout");
		super.printOptions(out);
		}

	
	@Override
	public int doWork(String[] args)
		{
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+"o:"))!=-1)
			{
			switch(c)
				{
				case 'o': this.setOutputFile(opt.getOptArg());break;
				default:
					{
					switch(handleOtherOptions(c, opt,args))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		return mainWork(opt.getOptInd(), args);
		}

	
	public static void main(String[] args)
		{
		new Vcf2Xml().instanceMain(args);
		}
}
