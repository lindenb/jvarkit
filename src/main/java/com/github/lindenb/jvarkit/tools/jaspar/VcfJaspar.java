/**
 * 
 */
package com.github.lindenb.jvarkit.tools.jaspar;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.util.CloserUtil;

import org.broad.tribble.readers.LineIterator;
import org.broad.tribble.readers.LineIteratorImpl;
import org.broad.tribble.readers.LineReader;
import org.broad.tribble.readers.LineReaderUtil;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.variant.vcf.VCFHeader;
import org.broadinstitute.variant.vcf.VCFHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineCount;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.AbstractVCFFilter2;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfJaspar extends AbstractVCFFilter2
	{
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private List<Matrix> jasparDb=new ArrayList<Matrix>();
	private double fraction_of_max=0.95;

	private VcfJaspar() {
		}
	

	@Override
	public String getProgramDescription() {
		return "Finds JASPAR profiles in VCF";
		}
	@Override
	protected String getOnlineDocUrl() {
		return "https://github.com/lindenb/jvarkit/wiki/VcfJaspar";
		}
	
	
	@Override
	protected void doWork(VcfIterator in, VariantContextWriter out)
			throws IOException
		{
		final String ATT="JASPAR";
		GenomicSequence genomicSequence=null;
		VCFHeader header=in.getHeader();
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"CmdLine",String.valueOf(getProgramCommandLine())));
		header.addMetaDataLine(new VCFHeaderLine(getClass().getSimpleName()+"Version",String.valueOf(getVersion())));

		header.addMetaDataLine(new VCFInfoHeaderLine(ATT,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String, 
				"Jaspar pattern overlapping: Format: (Name|length|Score/1000|pos|strand)"
				));
		out.writeHeader(header);
		while(in.hasNext())
			{
			VariantContext var=in.next();

			if(genomicSequence==null || !genomicSequence.getChrom().equals(var.getChr()))
				{
				info("Loading sequence "+var.getChr());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,var.getChr());
				}
			
			Set<String> hits=new HashSet<String>();
		
			for(Matrix matrix:this.jasparDb)
				{
					int start0=Math.max(0, var.getStart() - matrix.length());
					for(int y=start0;y<var.getStart() && y+matrix.length() <= genomicSequence.length();++y)
						{
						
						CharSequence forward=new SubSequence(genomicSequence,y,y+matrix.length());
						CharSequence revcomp=new RevCompCharSequence(forward);
						
						//run each strand
						for(int strand=0;strand<2;++strand)
							{
							double score= matrix.score(strand==0?forward:revcomp);
							if(score<=0) continue;
							
							if(score>= matrix.max()*this.fraction_of_max)
								{
								StringBuilder b=new StringBuilder("(");
								b.append(matrix.getName().replaceAll("[ \t;=]+", "/"));
								b.append("|");
								b.append(matrix.length());
								b.append("|");
								b.append((int)(1000.0*(score/matrix.max())));
								b.append("|");
								b.append(y+1);
								b.append("|");
								b.append(strand==0?'+':'-');
								b.append(")");
								hits.add(b.toString());
								break;
								}
							}
						
					}
				}
			if(hits.isEmpty())
				{
				out.add(var);
				continue;
				}
			VariantContextBuilder vcb=new VariantContextBuilder(var);
			vcb.attribute(ATT, hits.toArray(new String[hits.size()]));
			out.add(vcb.make());
			}
		}
	@Override
	public void printOptions(PrintStream out)
		{
		out.println(" -J (uri) jaspar PFM uri. required. example: http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt ");
		out.println(" -R (fasta) path to reference sequence indexed with picard. Required.");
		out.println(" -f (0<ratio<1) fraction of best score. default:"+this.fraction_of_max);
		super.printOptions(out);
		}
	
	@Override
	public int doWork(String[] args)
		{
		File fasta=null;
		String jasparUri=null;
		com.github.lindenb.jvarkit.util.cli.GetOpt opt=new com.github.lindenb.jvarkit.util.cli.GetOpt();
		int c;
		while((c=opt.getopt(args,getGetOptDefault()+ "J:R:f:"))!=-1)
			{
			switch(c)
				{
				case 'J': jasparUri=opt.getOptArg(); break;
				case 'R': fasta=new File(opt.getOptArg()); break;
				case 'f': this.fraction_of_max=Double.parseDouble(opt.getOptArg()); break;
				case ':': System.err.println("Missing argument for option -"+opt.getOptOpt());return -1;
				default:
					{
					switch(handleOtherOptions(c, opt))
						{
						case EXIT_FAILURE: return -1;
						case EXIT_SUCCESS: return 0;
						default:break;
						}
					}
				}
			}
		if(jasparUri==null)
			{
			error("Undefined jaspar-uri");
			return -1;
			}

		
	
		
		
		
		if(fasta==null)
			{
			super.error("Undefined fasta sequence");
			return -1;
			}
		try
			{
			info("Reading JASPAR: "+jasparUri);
			LineReader lr=LineReaderUtil.fromBufferedStream(IOUtils.openURIForReading(jasparUri));
			LineIterator liter=new LineIteratorImpl(lr);
			Iterator<Matrix> miter=Matrix.iterator(liter);
			while(miter.hasNext())
				{
				Matrix matrix = miter.next();
				this.jasparDb.add(matrix.convertToPWM());
				}
			lr.close();
			info("JASPAR size: "+this.jasparDb.size());

			if(jasparDb.isEmpty())
				{
				warning("JASPAR IS EMPTY");
				}
			
			info("opening "+fasta);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(fasta);
			return doWork(opt.getOptInd(), args);
			}
		catch(Exception err)
			{
			error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			}
		}

	public static void main(String[] args)
		{
		new VcfJaspar().instanceMainWithExit(args);
		}
	}
