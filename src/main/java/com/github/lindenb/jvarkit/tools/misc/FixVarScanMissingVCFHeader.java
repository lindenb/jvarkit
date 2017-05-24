package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.BufferedLineReader;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;

/**
BEGIN_DOC


END_DOC
*/
@Program(name="fixvarscanmissingheader ",description="Fix the sample name in the #CHROM header and fix VCF header missingh in varscan2. Sometimes, but not always, it happends. See  http://seqanswers.com/forums/showthread.php?t=33235 ")
public class FixVarScanMissingVCFHeader extends Launcher {
	private static Log LOG=Log.getInstance(FixVarScanMissingVCFHeader.class);


	
    @Parameter(names="-S",description="add this sample")
    public List<String> SAMPLES=new ArrayList<String>();
    
    private void header()
    	{
    	stdout().print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if(!this.SAMPLES.isEmpty())
			{
			stdout().print("\tFORMAT");
			for(String S:this.SAMPLES)
				{
				stdout().print("\t");
				stdout().print(S);
				}
			}
		stdout().println();
		}
    @Override
    public int doWork(List<String> args) {
		LOG.info("reading from stdin ... ");
		try
			{
			
			int c=stdin().read();
			if(c==-1 || c!='#')
				{ 
				LOG.info("VCF header missing (Sample are "+this.SAMPLES+") . Fixing.");
				stdout().println("##fileformat=VCFv4.1");
				stdout().println("##source=VarScan2");
				stdout().println("##varscan2header=missing");
				stdout().println("##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= 15\">");
				stdout().println("##INFO=<ID=WT,Number=1,Type=Integer,Description=\"Number of samples called reference (wild-type)\">");
				stdout().println("##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of samples called heterozygous-variant\">");
				stdout().println("##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Number of samples called homozygous-variant\">");
				stdout().println("##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">");
				stdout().println("##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">");
				stdout().println("##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">");
				stdout().println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
				stdout().println("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
				stdout().println("##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">");
				stdout().println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= 15\">");
				stdout().println("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">");
				stdout().println("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">");
				stdout().println("##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">");
				stdout().println("##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">");
				stdout().println("##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">");
				stdout().println("##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">");
				stdout().println("##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">");
				stdout().println("##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">");
				stdout().println("##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">");
				stdout().println("##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">");
				header();
				IOUtils.copyTo(System.in, stdout());
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
						stdout().println("##varscan2samples=replace");
						header();//problem with varscan:it doesn't know the sample names :-) !!
						}
					else
						{
						stdout().println(line);
						}
					}
				r.close();
				}
			else
				{
				throw new IOException("BAD varscan input: first letter is ascii("+c+")");
				}
			
			
			}
		catch(IOException err)
			{
			LOG.error(err,"Boum");
			return -1;
			}
		return 0;
		}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new FixVarScanMissingVCFHeader().instanceMainWithExit(args);

	}

}
