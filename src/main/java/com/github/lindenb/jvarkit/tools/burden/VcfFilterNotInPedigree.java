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
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

import com.github.lindenb.jvarkit.util.Pedigree;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

/**
 * VcfFilterNotInPedigree
 *
 */
public class VcfFilterNotInPedigree
	extends AbstractVcfFilterNotInPedigree
	{
	private static final org.slf4j.Logger LOG = com.github.lindenb.jvarkit.util.log.Logging.getLog(VcfFilterNotInPedigree.class);
	
	
	public VcfFilterNotInPedigree()
		{
		}
	 
	
	/* public for knime */
	@Override
	public Collection<Throwable> doVcfToVcf(
			final String inputName,
			final VcfIterator in,
			final VariantContextWriter out
			) throws IOException {
		final int IGNORE_SINGLETON=-1;
		final VCFHeader header = in.getHeader();
		
		try {
			final Pedigree pedigree = Pedigree.readPedigree(header);
			if(pedigree.isEmpty()) {
				return wrapException("No pedigree found in header "+inputName+". use VcfInjectPedigree to add it");
				}
			if(!pedigree.verifyPersonsHaveUniqueNames()) {
				return wrapException("I can't use this pedigree in VCF because two samples have the same ID in  "+inputName);
			}

			final Set<String> samplesNames= new HashSet<>(header.getSampleNamesInOrder());
			final Set<Pedigree.Person> individuals = new HashSet<>(pedigree.getPersons());
			final Iterator<Pedigree.Person> iter= individuals.iterator();
			while(iter.hasNext())
			{
				final Pedigree.Person person = iter.next();
				if(!(samplesNames.contains(person.getId()))) {
					LOG.warn("Ignoring "+person+" because not in VCF header or status is unknown");
					iter.remove();
				}
			}
			
			final VCFFilterHeaderLine filter = new VCFFilterHeaderLine(
					super.filterName,
					"Will be set for variant where the only genotypes non-homref are NOT in the pedigree "
					);
			final VCFFilterHeaderLine singletonFilter = new VCFFilterHeaderLine(
					super.singletonfilterName,
					"The ALT allele is found in less or equals than "+super.singleton+" individuals in the cases/controls"
					);
			
			final VCFHeader h2= new VCFHeader(header);
			h2.addMetaDataLine(filter);
			if(super.singleton!=IGNORE_SINGLETON) {
				h2.addMetaDataLine(singletonFilter);
			}
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary());
			out.writeHeader(h2);
			while(in.hasNext() &&  !out.checkError())
				{
				final VariantContext ctx = progess.watch(in.next());
				final boolean in_pedigree=individuals.stream().
						map(P->ctx.getGenotype(P.getId())).
						anyMatch(g->(!(g==null ||
								!g.isCalled() ||
								!g.isAvailable() ||
								g.isNoCall() ||
								g.isHomRef())))
								;
				
								
				if(!in_pedigree) {
					if(super.dicardVariant) continue;
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					vcb.filter(filter.getID());
					out.add(vcb.make());
					}
				else
					{
					boolean is_singleton;
					if(super.singleton!=IGNORE_SINGLETON) {
						is_singleton = true;
						for(final Allele alt:ctx.getAlternateAlleles()) {
							if( individuals.stream().
									map(P->ctx.getGenotype(P.getId())).
									filter(g->g.isCalled()&& g.getAlleles().contains(alt)).
									count() > super.singleton) 
								{
								is_singleton =false;
								break;
								}
							}
						}
					else
						{
						is_singleton=false;
						}
					if(is_singleton) {
						if(super.dicardVariant) continue;
						final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
						vcb.filter(singletonFilter.getID());
						out.add(vcb.make());
						}
					else
						{
						out.add(ctx);
						}
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
		new VcfFilterNotInPedigree().instanceMainWithExit(args);
		}
	}
