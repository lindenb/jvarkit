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
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * injects pedigree in VCF header
 *
 */
public class VcfInjectPedigree
	extends AbstractVcfInjectPedigree
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfInjectPedigree.class);
	
	
	public VcfInjectPedigree()
		{
		}
	 
	@Override
	public Collection<Throwable> initializeKnime() {
		if(super.pedigreeFile==null || !super.pedigreeFile.exists()) {
			return wrapException("Undefined Pedigree file option -"+OPTION_PEDIGREEFILE);
			}
		return super.initializeKnime();
	 	}
	
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out
			) throws IOException {
		final VCFHeader header = in.getHeader();
		final List<String> sampleNames = Collections.unmodifiableList(header.getSampleNamesInOrder());

		
		try {
			LOG.info("reading "+super.pedigreeFile);
			final Pedigree pedigree = Pedigree.readPedigree(super.pedigreeFile);
			
			if(!super.ignorePedigreeValidation)
				{
				pedigree.validate();
				}
			final Set<String> personsIds= new HashSet<>();
			for(final Pedigree.Person p:pedigree.getPersons()) {
				if(personsIds.contains(p.getId())) {
					return wrapException(
							"Pedigree contains two individual with ID :" +p.getId()
							+" This is not a pedigree error, but there is no mechanism to find them in the VCF header.");
					
				}
				personsIds.add(p.getId());
			}
			
			final Set<VCFHeaderLine> metaData = new HashSet<>(header.getMetaDataInInputOrder());
			
			/* remove previous pedigree */
			final Iterator<VCFHeaderLine> iter = metaData.iterator();
			while(iter.hasNext()) {
				final VCFHeaderLine hl = iter.next();
				if(Pedigree.VcfHeaderKey.equals( hl.getKey())) {
					if(super.cleanPreviousPedigree) {
					LOG.info("removing "+hl);
					iter.remove();
					}
				}
			}
			
			/* check already existing individual in vcf header */
			final Pedigree previousPedigree = Pedigree.readPedigree(metaData);
			for(final Pedigree.Family fam:pedigree.getFamilies()) {
				final Pedigree.Family prevFam = previousPedigree.getFamilyById(fam.getId());
				if(prevFam==null) continue;
				for(final Pedigree.Person person:fam.getIndividuals()) {
					final Pedigree.Person prevPers = prevFam.getPersonById(person.getId());
					if(prevPers!=null) {
						return wrapException("VCFHeader from "+inputName+" already contains  sample "+person+
								". Use option -"+OPTION_CLEANPREVIOUSPEDIGREE+" to removed it.");
					}
				}
			}
			
			/* check missing in header */
			for(final Pedigree.Person p: pedigree.getPersons()) {
				if(!sampleNames.contains(p.getId())) {
					if(!super.ignoreMissingInHeader) 
						{
						return wrapException("Sample "+p+" is declared in pedigree is missing in VCF header."+
								"Use option -"+OPTION_IGNOREMISSINGINHEADER+" to disable this error.");
						}
					else
						{
						LOG.warn("Sample "+p+" in pedigree is missing in VCF header.");
						}
					}
				}
			
			/* check missing in pedigree */
			for(final String sampleName: sampleNames) {
				if(!personsIds.contains(sampleName)) {
					if(!super.ignoreMissingInPedigree) 
					{
					return wrapException("VCF Sample "+sampleName+"  is missing in pedigree."+
							"Use option -"+OPTION_IGNOREMISSINGINPEDIGREE+" to disable this error.");
					}
				else
					{
					LOG.warn("Sample "+sampleName+" is present in VCF header is missing in Pedigree.");
					}
				}
			}
			
			
			metaData.addAll(pedigree.toVCFHeaderLines());
			addMetaData(metaData);
			
		
			final VCFHeader h2= new VCFHeader(metaData, sampleNames);
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				out.add(progess.watch(in.next()));
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
		new VcfInjectPedigree().instanceMainWithExit(args);
		}
	}
