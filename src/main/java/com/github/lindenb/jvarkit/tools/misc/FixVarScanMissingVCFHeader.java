package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.Collection;
import com.github.lindenb.jvarkit.util.command.Command;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.BufferedLineReader;

import com.github.lindenb.jvarkit.io.IOUtils;

public class FixVarScanMissingVCFHeader extends AbstractFixVarScanMissingVCFHeader
	{
	private static final org.apache.commons.logging.Log LOG = org.apache.commons.logging.LogFactory.getLog(FixVarScanMissingVCFHeader.class);

	@Override
	public  Command createCommand() {
			return new MyCommand();
		}
		 
	public  static class MyCommand extends AbstractFixVarScanMissingVCFHeader.AbstractFixVarScanMissingVCFHeaderCommand
	{
		    
	    private void header(PrintStream out)
	    	{
	    	out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			if(!this.SAMPLES.isEmpty())
				{
				out.print("\tFORMAT");
				for(String S:this.SAMPLES)
					{
					out.print("\t");
					out.print(S);
					}
				}
			out.println();
			}
    
	private void doWork(InputStream in,PrintStream out) throws IOException
		{
			int c= in.read();
			if(c==-1 || c!='#')
				{ 
				LOG.info("VCF header missing (Sample are "+this.SAMPLES+") . Fixing.");
				out.println("##fileformat=VCFv4.1");
				out.println("##source=VarScan2");
				out.println("##varscan2header=missing");
				out.println("##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= 15\">");
				out.println("##INFO=<ID=WT,Number=1,Type=Integer,Description=\"Number of samples called reference (wild-type)\">");
				out.println("##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of samples called heterozygous-variant\">");
				out.println("##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Number of samples called homozygous-variant\">");
				out.println("##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">");
				out.println("##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">");
				out.println("##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">");
				out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
				out.println("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
				out.println("##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">");
				out.println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= 15\">");
				out.println("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">");
				out.println("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">");
				out.println("##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">");
				out.println("##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">");
				out.println("##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">");
				out.println("##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">");
				out.println("##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">");
				out.println("##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">");
				out.println("##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">");
				out.println("##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">");
				header(out);
				IOUtils.copyTo(System.in, out);
				}
			else if(c=='#')
				{
				String line;
				BufferedLineReader r=new BufferedLineReader(System.in);
				while((line=r.readLine())!=null)
					{
					if(c!=-1)
						{
						line="#"+line;
						c=-1;
						}
					if(line.startsWith("#CHROM\t"))
						{
						out.println("##varscan2samples=replace");
						header(out);//problem with varscan:it doesn't know the sample names :-) !!
						}
					else
						{
						out.println(line);
						}
					}
				r.close();
				}
			else
				{
				throw new IOException("BAD varscan input: first letter is ascii("+c+")");
				}
			
		
		}

	@Override
	protected Collection<Throwable> call(String inputName) throws Exception
		{
		PrintStream out = null; 
		InputStream in = null;
		try {
			in = (inputName==null?stdin():IOUtils.openURIForReading(inputName));
			out = openFileOrStdoutAsPrintStream();
			doWork(in,out);
			out.flush();
			return RETURN_OK;
		} catch (Exception e) {
			return wrapException(e);
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FixVarScanMissingVCFHeader().instanceMainWithExit(args);

	}

}
