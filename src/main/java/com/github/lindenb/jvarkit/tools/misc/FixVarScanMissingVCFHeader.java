package com.github.lindenb.jvarkit.tools.misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.github.lindenb.jvarkit.util.picard.cmdline.Option;
import com.github.lindenb.jvarkit.util.picard.cmdline.Usage;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.BufferedLineReader;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.picard.AbstractCommandLineProgram;

public class FixVarScanMissingVCFHeader extends AbstractCommandLineProgram {
	private static Log LOG=Log.getInstance(FixVarScanMissingVCFHeader.class);

    @Usage(programVersion="1.0")
    public String USAGE = getStandardUsagePreamble() + " Fix the sample name in the #CHROM header and fix VCF header missingh in varscan2. Sometimes, but not always, it happends. See  http://seqanswers.com/forums/showthread.php?t=33235 ";

	
    @Option(shortName="S",doc="add this sample",minElements=0)
    public List<String> SAMPLES=new ArrayList<String>();
    
    private void header()
    	{
    	System.out.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
		if(!this.SAMPLES.isEmpty())
			{
			System.out.print("\tFORMAT");
			for(String S:this.SAMPLES)
				{
				System.out.print("\t");
				System.out.print(S);
				}
			}
		System.out.println();
		}
    
	@Override
	protected int doWork()
		{
		LOG.info("reading from stdin ... ");
		try
			{
			
			int c=System.in.read();
			if(c==-1 || c!='#')
				{ 
				LOG.info("VCF header missing (Sample are "+this.SAMPLES+") . Fixing.");
				System.out.println("##fileformat=VCFv4.1");
				System.out.println("##source=VarScan2");
				System.out.println("##varscan2header=missing");
				System.out.println("##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= 15\">");
				System.out.println("##INFO=<ID=WT,Number=1,Type=Integer,Description=\"Number of samples called reference (wild-type)\">");
				System.out.println("##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of samples called heterozygous-variant\">");
				System.out.println("##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Number of samples called homozygous-variant\">");
				System.out.println("##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">");
				System.out.println("##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">");
				System.out.println("##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">");
				System.out.println("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
				System.out.println("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
				System.out.println("##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">");
				System.out.println("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= 15\">");
				System.out.println("##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">");
				System.out.println("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">");
				System.out.println("##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">");
				System.out.println("##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">");
				System.out.println("##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">");
				System.out.println("##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">");
				System.out.println("##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">");
				System.out.println("##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">");
				System.out.println("##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">");
				System.out.println("##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">");
				header();
				IOUtils.copyTo(System.in, System.out);
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
						System.out.println("##varscan2samples=replace");
						header();//problem with varscan:it doesn't know the sample names :-) !!
						}
					else
						{
						System.out.println(line);
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
