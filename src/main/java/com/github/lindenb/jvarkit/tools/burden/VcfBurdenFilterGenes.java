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

*/
package com.github.lindenb.jvarkit.tools.burden;

import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParser;
import com.github.lindenb.jvarkit.util.vcf.predictions.VepPredictionParserFactory;

/**
 * Filter by genes
 *
 */
public class VcfBurdenFilterGenes
	extends AbstractVcfBurdenFilterGenes
	{
	private final Set<String> geneNames= new HashSet<>();
	
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfBurdenFilterGenes.class);
	
	
	public VcfBurdenFilterGenes()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.geneFile==null || !super.geneFile.exists()) {
		return wrapException("Undefined gene file option -"+OPTION_GENEFILE);
		}
		
		try {
			this.geneNames.clear();
			this.geneNames.addAll(Files.readAllLines(super.geneFile.toPath()));
		} catch (IOException e) {
			return wrapException(e);
		}
		geneNames.remove(".");
		geneNames.remove("");
		LOG.info("number of genes : "+geneNames.size());
		return super.initializeKnime();
	 	}

	@Override
	public void disposeKnime() {
		this.geneNames.clear();
		super.disposeKnime();
	}
	
	
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out
			) throws IOException {
		final VCFHeader header=in.getHeader();		
		try {
			final VCFHeader h2=addMetaData(new VCFHeader(header));

			final VCFFilterHeaderLine filterControlsHeader;
			if(!super.filterTag.trim().isEmpty())
				{
				filterControlsHeader = new VCFFilterHeaderLine(
						super.filterTag.trim(),
						"Genes in list "+this.filterTag
						);
				h2.addMetaDataLine(filterControlsHeader);
				}
			else
				{
				filterControlsHeader = null;
				}
			final List<String> lookColumns=Arrays.asList(
					"CCDS",
					"Feature",
					"ENSP",
					"Gene",
					"HGNC",
					"HGNC_ID",					
					"SYMBOL",					
					"RefSeq"					
					);
			final VepPredictionParser vepParser = new VepPredictionParserFactory(header).get();
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.watch(in.next());
				boolean keep=false;
				final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
				
				
				final List<String> csqList= ctx.getAttributeAsList(vepParser.getTag()).stream().
						map(O->String.class.cast(O)).collect(Collectors.toList());
				
				//not just set FILTER ?
				if(filterControlsHeader==null) {
					vcb.rmAttribute(vepParser.getTag());
				}
				
				final List<String> newCsqList=new ArrayList<>();
				for(final String predStr: csqList) {
					final VepPredictionParser.VepPrediction pred = vepParser.parseOnePrediction(ctx,predStr);
					for(final String col:lookColumns) {
						final String token = pred.getByCol(col);
						if(token!=null && !token.isEmpty() && this.geneNames.contains(token))
							{
							newCsqList.add(predStr);
							keep=true;
							break;
							}
						}
					}
				//not just set FILTER ?
				if(filterControlsHeader==null) {
					vcb.attribute(vepParser.getTag(),newCsqList);
				}
				
				
				if(filterControlsHeader!=null)
					{
					if(!keep) vcb.filter(filterControlsHeader.getID());
					out.add(ctx);
					}
				else
					{
					if(keep) out.add(vcb.make());
					}
				}
			progess.finish();
			return RETURN_OK;
			} catch(Exception err) {
				return wrapException(err);
			} finally {
				CloserUtil.close(in);
			}
		}
	
	@Override
	protected Collection<Throwable> call(String inputName) throws Exception {
		return doVcfToVcf(inputName);
		}
	 	
	
	public static void main(String[] args)
		{
		new VcfBurdenFilterGenes().instanceMainWithExit(args);
		}
	}
