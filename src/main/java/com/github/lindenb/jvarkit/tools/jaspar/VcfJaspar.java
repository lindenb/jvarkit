/**
 * 
 */
package com.github.lindenb.jvarkit.tools.jaspar;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.CloserUtil;

import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.SubSequence;
import com.github.lindenb.jvarkit.util.bio.RevCompCharSequence;
import com.github.lindenb.jvarkit.util.picard.GenomicSequence;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;


public class VcfJaspar extends AbstractVcfJaspar
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfJaspar.class);
	private IndexedFastaSequenceFile indexedFastaSequenceFile=null;
	private List<Matrix> jasparDb=new ArrayList<Matrix>();
	

	public VcfJaspar() {
		}
	

	@Override
	protected Collection<Throwable> doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out)
			throws IOException {
		final String ATT="JASPAR";
		GenomicSequence genomicSequence=null;
		final VCFHeader header=new VCFHeader(in.getHeader());
		addMetaData(header);

		header.addMetaDataLine(new VCFInfoHeaderLine(ATT,
				VCFHeaderLineCount.UNBOUNDED,
				VCFHeaderLineType.String, 
				"Jaspar pattern overlapping: Format: (Name|length|Score/1000|pos|strand)"
				));
		out.writeHeader(header);
		while(in.hasNext())
			{
			VariantContext var=in.next();

			if(genomicSequence==null || !genomicSequence.getChrom().equals(var.getContig()))
				{
				LOG.info("Loading sequence "+var.getContig());
				genomicSequence=new GenomicSequence(this.indexedFastaSequenceFile,var.getContig());
				}
			
			final Set<String> hits=new HashSet<String>();
		
			for(final Matrix matrix:this.jasparDb)
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
			final VariantContextBuilder vcb=new VariantContextBuilder(var);
			vcb.attribute(ATT, hits.toArray(new String[hits.size()]));
			out.add(vcb.make());
			}
		return RETURN_OK;
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		if(super.jasparUri==null)
			{
			return wrapException("Undefined jaspar-uri");
			}		
		
		if(super.fasta==null)
			{
			return wrapException("Undefined fasta sequence");
			}
		try
			{
			LOG.info("Reading JASPAR: "+jasparUri);
			LineIterator liter = IOUtils.openURIForLineIterator(this.jasparUri);
			Iterator<Matrix> miter=Matrix.iterator(liter);
			while(miter.hasNext())
				{
				Matrix matrix = miter.next();
				this.jasparDb.add(matrix.convertToPWM());
				}
			CloserUtil.close(liter);
			LOG.info("JASPAR size: "+this.jasparDb.size());

			if(jasparDb.isEmpty())
				{
				LOG.warn("JASPAR IS EMPTY");
				}
			
			LOG.info("opening "+fasta);
			this.indexedFastaSequenceFile=new IndexedFastaSequenceFile(fasta);
			return doVcfToVcf(inputName);
			}
		catch(Exception err)
			{
			return wrapException(err);
			}
		finally
			{
			CloserUtil.close(this.indexedFastaSequenceFile);
			this.indexedFastaSequenceFile=null;
			}
		}

	public static void main(String[] args)
		{
		new VcfJaspar().instanceMainWithExit(args);
		}
	}
