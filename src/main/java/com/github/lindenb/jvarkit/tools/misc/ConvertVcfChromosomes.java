/*
The MIT License (MIT)

Copyright (c) 2016 Pierre Lindenbaum

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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

public class ConvertVcfChromosomes extends AbstractConvertVcfChromosomes {
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(ConvertVcfChromosomes.class);
	private Map<String,String> customMapping=new HashMap<String,String>();
	private Set<String> unmappedChromosomes=new HashSet<String>();

	public ConvertVcfChromosomes()
		{
		}
	
	private String convertName(final String chrom)throws IOException
		{
		if(chrom==null) throw new NullPointerException();
		String newname=customMapping.get(chrom);
		if(newname==null)
			{
			if(!unmappedChromosomes.contains(chrom))
				{
				LOG.warn("unmapped chromosome "+chrom);
				unmappedChromosomes.add(chrom);
				}
			if(ignore_if_no_mapping) return null;
			
			if(use_original_chrom_name_if_no_mapping)
				{	
				return chrom;
				}
			throw new IOException("No mapping found to convert name of chromosome \""+chrom+"\"");
			}
		return newname;
		}
	
	/** public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final  VcfIterator in,
			final  VariantContextWriter out)
			throws IOException {
		final VCFHeader header1=in.getHeader();
	
		final  Set<VCFHeaderLine> meta2=new LinkedHashSet<VCFHeaderLine>();
		for(final  VCFHeaderLine L:header1.getMetaDataInInputOrder())
			{
			if(!L.getKey().equals(VCFHeader.CONTIG_KEY))
				{
				meta2.add(L);
				}
			}
		final VCFHeader header2=new VCFHeader(meta2,header1.getSampleNamesInOrder());

		if(header1.getSequenceDictionary()!=null)
			{
			final List<SAMSequenceRecord> ssrs=new ArrayList<SAMSequenceRecord>();

			for(int i=0;i< header1.getSequenceDictionary().size();++i)
				{
				SAMSequenceRecord ssr=header1.getSequenceDictionary().getSequence(i);
				String newName=convertName(ssr.getSequenceName());
				if(newName==null)
					{
					//skip unknown chromosomes
					continue;
					}
				ssr=new SAMSequenceRecord(newName, ssr.getSequenceLength());
				ssrs.add(ssr);
				}
			header2.setSequenceDictionary(new SAMSequenceDictionary(ssrs));
			}
		super.addMetaData(header2);
		out.writeHeader(header2);
		
		while(in.hasNext())
			{
			final VariantContext ctx=in.next();
			final String newName=convertName(ctx.getContig());
			if(newName==null)
				{
				//skip unknown chromosomes
				continue;
				}
			final VariantContextBuilder vcb= super.getVariantContextBuilderFactory().newVariantContextBuilder(ctx);
			vcb.chr(newName);
			out.add(vcb.make());
			}
		
		if(!unmappedChromosomes.isEmpty())
			{
			LOG.warn("Unmapped chromosomes: "+unmappedChromosomes);
			}
		return RETURN_OK;
		}
	
	@Override
	public Collection<Throwable> initializeKnime() {
		if( super.mappingFile==null)
			{
			return wrapException("undefined mapping file");
			}
		IOUtil.assertFileIsReadable(super.mappingFile);
		this.customMapping.clear();
		BufferedReader in=null;
		try
			{
			LOG.info("Loading custom mapping "+super.mappingFile);
			in=IOUtils.openFileForBufferedReading(super.mappingFile);
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.isEmpty() || line.startsWith("#")) continue;
				String tokens[]=line.split("[\t]");
				if(tokens.length!=2
						|| tokens[0].trim().isEmpty()
						|| tokens[1].trim().isEmpty()
						) {
					in.close();
					throw new IOException("Bad mapping line: \""+line+"\"");
					}
				tokens[0]=tokens[0].trim();
				tokens[1]=tokens[1].trim();
				if(this.customMapping.containsKey(tokens[0]))
					{
					in.close();
					throw new IOException("Mapping defined twice for: \""+tokens[0]+"\"");
					}
				this.customMapping.put(tokens[0], tokens[1]);
				}
			}
		catch(IOException err) {
			return wrapException(err);
		}
		finally
			{
			CloserUtil.close(in);
			}
		return super.initializeKnime();
	}
	
	@Override
	public void disposeKnime() {
		this.customMapping.clear();
		super.disposeKnime();
		}
	
	@Override
	protected Collection<Throwable> call(final String inputName) throws Exception {
		return super.doVcfToVcf(inputName);
		}
	

	/**
	 * @param args
	 */
	public static void main(String[] args)
		{
		new ConvertVcfChromosomes().instanceMainWithExit(args);
		}
	}
