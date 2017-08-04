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

import java.io.File;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
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
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC


END_DOC
 */
@Program(name="vcffilternotinpedigree",
	description="Adds a FILTER NotInPedigree if the only not(homref) genotypes are not in a pedigree",
	keywords={"burden","vcf","pedigree"}
	)
public class VcfFilterNotInPedigree
	extends Launcher
	{

	private static final Logger LOG = Logger.build(VcfFilterNotInPedigree.class).make();


	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;


	@Parameter(names={"-f","--filter"},description="FILTER name. Will be set for variant where the only genotypes non-homref are NOT in the pedigree")
	private String filterName = "NoGenotypeInPedigree";

	@Parameter(names={"-r","--remove"},description="remove the variant instead of setting the FILTER")
	private boolean dicardVariant = false;

	@Parameter(names={"-s","--singleton"},description="Variant is flagged/FILTERed as SingletonAlt if the ALT if found in less or equal times 'singleton-times' in the genotypes. -1:ignore")
	private int singleton = 1 ;

	@Parameter(names={"-sf","--sfilter"},description="FILTER name for option singleton")
	private String singletonfilterName = "SingletonAlt";
	
	
	public VcfFilterNotInPedigree()
		{
		}
	 
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
	final int IGNORE_SINGLETON=-1;
		final VCFHeader header = in.getHeader();
		
		try {
			final Pedigree pedigree = Pedigree.newParser().parse(header);
			if(pedigree.isEmpty()) {
				LOG.error("No pedigree found in header "+inputName+". use VcfInjectPedigree to add it");
				return -1;
				}
			if(!pedigree.verifyPersonsHaveUniqueNames()) {
				LOG.error("I can't use this pedigree in VCF because two samples have the same ID in  "+inputName);
				return -1;
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
					this.filterName,
					"Will be set for variant where the only genotypes non-homref are NOT in the pedigree "
					);
			final VCFFilterHeaderLine singletonFilter = new VCFFilterHeaderLine(
					this.singletonfilterName,
					"The ALT allele is found in less or equals than "+this.singleton+" individuals in the cases/controls"
					);
			
			final VCFHeader h2= new VCFHeader(header);
			h2.addMetaDataLine(filter);
			if(this.singleton!=IGNORE_SINGLETON) {
				h2.addMetaDataLine(singletonFilter);
			}
			
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header).logger(LOG);
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
					if(this.dicardVariant) continue;
					final VariantContextBuilder vcb = new VariantContextBuilder(ctx);
					vcb.filter(filter.getID());
					out.add(vcb.make());
					}
				else
					{
					boolean is_singleton;
					if(this.singleton!=IGNORE_SINGLETON) {
						is_singleton = true;
						for(final Allele alt:ctx.getAlternateAlleles()) {
							if( individuals.stream().
									map(P->ctx.getGenotype(P.getId())).
									filter(g->g.isCalled()&& g.getAlleles().contains(alt)).
									count() > this.singleton) 
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
						if(this.dicardVariant) continue;
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
				LOG.error(err);
				return -1;
			} finally {
				CloserUtil.close(in);
			}
		}
	
	
	@Override
	public int doWork(final List<String> args) {
		return doVcfToVcf(args, outputFile);
	}
	
	
	public static void main(String[] args)
		{
		new VcfFilterNotInPedigree().instanceMainWithExit(args);
		}
	}
