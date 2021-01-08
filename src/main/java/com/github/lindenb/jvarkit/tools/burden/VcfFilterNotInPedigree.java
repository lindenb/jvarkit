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
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.jcommander.OnePassVcfLauncher;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.pedigree.Pedigree;
import com.github.lindenb.jvarkit.pedigree.PedigreeParser;
import com.github.lindenb.jvarkit.pedigree.Sample;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;

/**
BEGIN_DOC

## Example



END_DOC
 */
@Program(name="vcffilternotinpedigree",
	description="Adds a FILTER 'NotInPedigree' if the only not(homref) genotypes are not in a pedigree",
	keywords={"burden","vcf","pedigree"},
	modificationDate="20200713"
	)
public class VcfFilterNotInPedigree
	extends OnePassVcfLauncher
	{

	private static final Logger LOG = Logger.build(VcfFilterNotInPedigree.class).make();
	@Parameter(names={"-f","--filter"},description="FILTER name. Will be set for variant where the only genotypes non-homref are NOT in the pedigree")
	private String filterName = "NoGenotypeInPedigree";

	@Parameter(names={"-r","--remove"},description="remove the variant instead of setting the FILTER (hard filtering)")
	private boolean dicardVariant = false;

	@Parameter(names={"-s","--singleton"},description="Variant is flagged/FILTERed as SingletonAlt if the ALT if found in less or equal times 'singleton-times' in the genotypes. -1:ignore")
	private int singleton = 1 ;
	
	@Parameter(names={"-sf","--sfilter"},description="FILTER name for option --singleton")
	private String singletonfilterName = "SingletonAlt";
	
	@Parameter(names={"-p","--pedigree"},description=PedigreeParser.OPT_DESC,required=true)
	private Path pedigreeFile = null;

	@Parameter(names={"-gtf","--ignore-filtered-gt"},description="[20180406] Do not consider a *genotype* if it is FILTERED.")
	private boolean ignoreFilteredGT = false;

	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	
	@Override
	protected int doVcfToVcf(String inputName, VCFIterator iterin, VariantContextWriter out) {
		final int IGNORE_SINGLETON=-1;
		final Predicate<Genotype> acceptFilteredGenotype = G -> 
			G!=null && (this.ignoreFilteredGT == false ||  !G.isFiltered());
			final Pedigree pedigree;
		final VCFHeader header = iterin.getHeader();
		
			
	
			try {
				pedigree = new PedigreeParser().parse(this.pedigreeFile);
			} catch (final IOException err) {
				LOG.error("Cannot read pedigree in file: "+this.pedigreeFile,err);
				return -1;
				}
			
			if(pedigree.isEmpty()) {
				throw new JvarkitException.PedigreeError("No pedigree found in header. use VcfInjectPedigree to add it");
				}
			pedigree.checkUniqIds();
			

			final Set<String> samplesNames= new HashSet<>(header.getSampleNamesInOrder());
			final Set<Sample> individuals = new HashSet<>(pedigree.getSamples());
			final Iterator<Sample> iter= individuals.iterator();
			while(iter.hasNext()) {
				final Sample person = iter.next();
				if(!(samplesNames.contains(person.getId()))) {
					LOG.warn("Ignoring "+person+" because not in VCF header.");
					iter.remove();
					}
				}
			
			final VCFFilterHeaderLine filter = new VCFFilterHeaderLine(
					this.filterName,
					"Will be set for variant where the only genotypes non-homref are NOT in the pedigree "
					);
			final VCFFilterHeaderLine singletonFilter = new VCFFilterHeaderLine(
					this.singletonfilterName,
					"The ALT allele is found in less or equals than "+ this.singleton+" individuals in the cases/controls"
					);
			
			final VCFHeader h2= new VCFHeader(header);
			JVarkitVersion.getInstance().addMetaData(this, h2);
			if(this.singleton!=IGNORE_SINGLETON) {
				h2.addMetaDataLine(singletonFilter);
				}					
			out.writeHeader(h2);
	
			while(iterin.hasNext()) {
				final VariantContext ctx = iterin.next();
				final boolean in_pedigree= individuals.stream().
						map(P->ctx.getGenotype(P.getId())).
						filter(acceptFilteredGenotype).
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
							if(individuals.stream().
									map(P->ctx.getGenotype(P.getId())).
									filter(acceptFilteredGenotype).
									filter(g->g!=null && g.isCalled()&& g.getAlleles().contains(alt)).
									count() >this.singleton) 
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
		return 0;
		}
	
	
	public static void main(final String[] args)
		{
		new VcfFilterNotInPedigree().instanceMainWithExit(args);
		}
	}
