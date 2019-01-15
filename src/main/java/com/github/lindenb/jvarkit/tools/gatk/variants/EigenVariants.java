/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum

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
* 2017 creation

*/
package com.github.lindenb.jvarkit.tools.gatk.variants;

import java.io.File;
import java.util.Map;

import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;

import com.github.lindenb.jvarkit.tools.vcfeigen.EigenInfoAnnotator;

import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

@DocumentedGATKFeature(
		summary="Variant Annotator for the data of https://xioniti01.u.hpc.mssm.edu/v1.1/ : Eigen makes use of a variety of functional annotations in both coding and noncoding regions (such as made available by the ENCODE and Roadmap Epigenomics projects), and combines them into one single measure of functional importance.",
		groupName = HelpConstants.DOCS_CAT_VARMANIP,
		extraDocs = {CommandLineGATK.class}
		)
public class EigenVariants extends AbstractVariantProcessor {
@Input(shortName="eigen",fullName="eigenDirectory",required=true,doc="The Eigen directory containing the tabix indexed files  *.tab.gz.")
protected File eigenDir=null;
@Argument(fullName="tabixPrefix",shortName="tabixPrefix",required=false,doc="Override prefix of tabix file in the eigen directory. Leave null for default.")
public String tabixPrefix = null;

private EigenInfoAnnotator annotator = null;

public EigenVariants() {
	}

@Override
public void initialize() {
	logger.info("opening eigen directory :"+eigenDir);	
	IOUtil.assertDirectoryIsReadable(eigenDir);
	this.annotator = new EigenInfoAnnotator(this.eigenDir);
	if( this.tabixPrefix!=null) {
		logger.info("changing eigen file prefix from  :"+this.annotator.getTabixPrefix()+" to "+this.tabixPrefix);	
		this.annotator.setTabixFilePrefix(this.tabixPrefix);
	}
	
	final VCFHeader header1= super.getVcfHeader();
	
	final  VCFHeader h2 = new VCFHeader(header1);
	for(final VCFInfoHeaderLine vihl: this.annotator.getInfoHeaderLines()) {
		if(h2.getInfoHeaderLine(vihl.getID())!=null) {
			throw new UserException.MalformedVCFHeader("VCF INFO "+vihl.getID()+" already defined in input VCF.");
		}
		h2.addMetaDataLine(vihl);
		}
	super.writer.writeHeader(h2);
	super.initialize();
	}

@Override
protected VariantContext mapVariant(
		final VariantContext ctx,
		final RefMetaDataTracker tracker,
		final ReferenceContext ref,
		final AlignmentContext context
		) {
	final Map<String,Object> m  = this.annotator.getAnnotations(ctx);
	if(m==null || m.isEmpty())
		{
		return ctx;
		}
	else
		{
		final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
		for(final String key: m.keySet()) {
			vcb.attribute(key, m.get(key));
			}
		return vcb.make();
		}
	}

@Override
public void onTraversalDone(Long result) {
	this.annotator.close();
	this.annotator=null;
	super.onTraversalDone(result);
	}
}
