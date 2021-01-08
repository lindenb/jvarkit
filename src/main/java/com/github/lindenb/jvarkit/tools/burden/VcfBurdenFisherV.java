/*
The MIT License (MIT)

Copyright (c) 2021 Pierre Lindenbaum

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
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamWriter;

import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.math.stats.FisherExactTest;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;

/**

BEGIN_DOC

Variant in that VCF should have one and only one ALT allele. Use https://github.com/lindenb/jvarkit/wiki/VcfMultiToOneAllele if needed.


END_DOC
*/

@Program(
		name="vcfburdenfisherv",
		description="Fisher Case / Controls per Variant (Vertical)",
		keywords={"vcf","burden","fisher"},
		creationDate="20160418",
		modificationDate="20200304"
		)
public class VcfBurdenFisherV
	extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfBurdenFisherV.class).make();
	public static final String VCF_HEADER_FISHER_VALUE="VCFBurdenFisherV";

	
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-p","--pedigree"},description= PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile=null;
	@Parameter(names={"-if","--ignore-filter"},description="accept variants having a FILTER column. Default is ignore variants with a FILTER column")
	private boolean acceptFiltered = false;
	@Parameter(names={"-table","--table"},description="Write statistics into that X-HTML file instead of inserting a header line the VCF (faster)")
	private Path tableOut = null;
	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();
	
	public VcfBurdenFisherV()
		{
		}
	 
	@Override
	protected int doVcfToVcf(
		final String inputName,
		final VCFIterator in,
		final VariantContextWriter out)
			{
			final boolean HAS_VARIANT = true;
			final boolean NO_VARIANT = false;
			Path tmVcfOut = null;
			VariantContextWriter tmpw = null;
			try {
				
			final VCFHeader header = in.getHeader();
			
			if(tableOut==null && header.getMetaDataLine(VCF_HEADER_FISHER_VALUE)!=null)
				{
				throw new JvarkitException.UserError(
						"VCF Header "+VCF_HEADER_FISHER_VALUE+" already specified in input");
				}
			
			
			
			
			final VCFHeader header2 = new VCFHeader(header);

						
			final Set<Sample> persons = new PedigreeParser().
					parse(this.pedigreeFile).
					getSamplesInVcfHeader(header).
					filter(S->S.isStatusSet()).
					collect(Collectors.toSet());
			
			if(persons.isEmpty()) {
				LOG.warn("No sample in pedigree + vcf header");
				return -1;
			}
			
			final Map<Sample,Boolean> indi2supervariant = new HashMap<>(persons.size());

			for(final Sample  person: persons) {
				indi2supervariant.put(person,NO_VARIANT);
				}
			
			if(this.tableOut!=null) {
				tmVcfOut  = null;
				}
			else if(this.outputFile==null) {
				tmVcfOut = Files.createTempFile("tmp.", FileExtensions.BCF);	
				} 
			else  {
				tmVcfOut  = Files.createTempFile(this.outputFile.getParent(),"tmp.", FileExtensions.BCF);
				}

			
			if(tmVcfOut!=null) {
				final VariantContextWriterBuilder vcwb=new VariantContextWriterBuilder();
				vcwb.setCreateMD5(false);
				vcwb.setReferenceDictionary(SequenceDictionaryUtils.extractRequired(header));
				vcwb.clearOptions();
				vcwb.setOutputPath(tmVcfOut);
				tmpw=vcwb.build();
				tmpw.writeHeader(header);
				}
			else
				{
				tmpw = null;
				}

			long count_variants = 0L;
			while(in.hasNext()) {
				final VariantContext ctx = in.next();
				
				if(tmpw!=null) {
					tmpw.add(ctx);
					}
				else
					{
					out.add(ctx);
					}
				
				if(ctx.isFiltered() && !this.acceptFiltered) continue;
				final int n_alts = ctx.getAlternateAlleles().size();
				
				if( n_alts == 0) {
					LOG.warn("ignoring variant without ALT allele.");
					continue;
				}
				
				count_variants++;
				
				if( n_alts > 1) {
					LOG.warn("variant with more than one ALT. "+ctx.getContig()+":"+ctx.getStart());
					}
				
				
				//loop over person in this pedigree
				for(final Sample person : indi2supervariant.keySet() ) {
					if(indi2supervariant.get(person)==HAS_VARIANT) continue;
					final Genotype g = ctx.getGenotype(person.getId());	
					if(g==null) {
						continue;//not in vcf header
						}
					if(g.isFiltered()) {
						LOG.warn("ignoring filtered genotype");
						continue;//not filter.
						}
					if(g.getAlleles().stream().anyMatch(A->A.isCalled() && A.isNonReference())) {
						indi2supervariant.put(person,HAS_VARIANT);
						}//end of allele
					}//en dof for[person]
			} //end of 
				
			
			int count_case_sv0 = 0;
			int count_ctrl_sv0 = 0;
			int count_case_sv1 = 0;
			int count_ctrl_sv1 = 0;
		
		
			for(final Sample person :indi2supervariant.keySet() ) {
				final boolean hasVariant = indi2supervariant.get(person);
				if(!hasVariant) {
					if(person.isAffected()) count_case_sv0++;
					else count_ctrl_sv0++;
				} else // AT_LEAST_ONE_VARIANT 
					{
					if(person.isAffected()) count_case_sv1++;
					else count_ctrl_sv1++;
					}
			}//end of person
		
		
		
		
			final FisherExactTest fisher = FisherExactTest.compute(
					count_case_sv0, count_case_sv1,
					count_ctrl_sv0, count_ctrl_sv1
					);
			
			
			if(this.tableOut==null) {
				header2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_FISHER_VALUE,
						String.valueOf(fisher.getAsDouble())));
				header2.addMetaDataLine(new VCFHeaderLine(VCF_HEADER_FISHER_VALUE+".count",
						String.join("|",
						"CASE_SV0="+count_case_sv0,
						"CASE_SV1="+count_case_sv1,
						"CTRL_SV0="+count_ctrl_sv0,
						"CTRL_SV1="+count_ctrl_sv1		
						)));
				try(VCFIterator in2 = new VCFIteratorBuilder().open(tmVcfOut)) {
					out.writeHeader(header2);
					while(in2.hasNext()) {
						out.add(in2.next());
						}
					}
				Files.deleteIfExists(tmVcfOut);
				tmpw=null;
				}
			else
				{
				/* save xml report */
				try(final OutputStream os=super.openPathOrStdoutAsPrintStream(this.tableOut)) {
					final XMLOutputFactory xof = XMLOutputFactory.newFactory();
					final XMLStreamWriter w = xof.createXMLStreamWriter(os, "UTF-8");
					w.writeStartDocument("UTF-8", "1.0");
					w.writeStartElement("html");
					w.writeStartElement("head");
					w.writeStartElement("title");
					w.writeCharacters(getClass().getSimpleName()+":"+inputName);
					w.writeEndElement();
					
					w.writeEmptyElement("meta");
					w.writeAttribute("name", "vcf");
					w.writeAttribute("content",String.valueOf(inputName));

					w.writeEmptyElement("meta");
					w.writeAttribute("name", "version");
					w.writeAttribute("content",JVarkitVersion.getInstance().getLabel());

					
					w.writeEndElement();//hread
					w.writeStartElement("body");
					w.writeStartElement("table");
					w.writeStartElement("caption");
					w.writeCharacters("Fisher: ");
					w.writeStartElement("span");
					w.writeAttribute("id", "fisher");
					w.writeCharacters(String.valueOf(fisher.getAsDouble()));
					w.writeEndElement();//span
					w.writeCharacters(" Variant(s): ");
					w.writeStartElement("span");
					w.writeAttribute("id", "variants");
					w.writeCharacters(String.valueOf(count_variants));
					w.writeEndElement();//span
					w.writeEndElement();//caption
					
					w.writeStartElement("tr");
					w.writeEmptyElement("th");
					w.writeStartElement("th");
					w.writeCharacters("With Rare");
					w.writeEndElement();//th
					w.writeStartElement("th");
					w.writeCharacters("No Rare");
					w.writeEndElement();//th
					w.writeEndElement();//tr
					
					w.writeStartElement("tr");
					w.writeStartElement("th");
					w.writeCharacters("Case");
					w.writeEndElement();//th
					w.writeStartElement("td");
					w.writeAttribute("id", "case1");
					w.writeCharacters(String.valueOf(count_case_sv1));
					w.writeEndElement();//td
					w.writeStartElement("td");
					w.writeAttribute("id", "case0");
					w.writeCharacters(String.valueOf(count_case_sv0));
					w.writeEndElement();//td
					w.writeEndElement();//tr

					w.writeStartElement("tr");
					w.writeStartElement("th");
					w.writeCharacters("Controls");
					w.writeEndElement();//th
					w.writeStartElement("td");
					w.writeAttribute("id", "ctrl1");
					w.writeCharacters(String.valueOf(count_ctrl_sv1));
					w.writeEndElement();//td
					w.writeStartElement("td");
					w.writeAttribute("id", "ctrl0");
					w.writeCharacters(String.valueOf(count_ctrl_sv0));
					w.writeEndElement();//td
					w.writeEndElement();//tr

					w.writeEndElement();//table
					w.writeEndElement();//body
					w.writeEndElement();//html
					w.writeEndDocument();
					w.close();
					os.flush();
					}
				}
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally {
			if(tmVcfOut!=null) try {Files.deleteIfExists(tmVcfOut);} catch(IOException err) {}
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		try 
			{
			return doVcfToVcfPath(args,this.writingVariantsDelegate,this.outputFile);
			}
		catch(final Exception err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args)
		{
		new VcfBurdenFisherV().instanceMainWithExit(args);
		}
	}
