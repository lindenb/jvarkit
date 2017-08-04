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

BEGIN_DOC


This tools reads a pedigree file and inject it in the VCF header  



```
$ java -jar dist/vcfinjectpedigree.jar \
	-imih -imip -p input.ped \
	input.vcf.gz > out.vcf

$ grep Sample out.vcf
(...)
##Sample=<Family=F1,ID=INDI1,Father=0,Mother=0,Sex=1,Status=1>
##Sample=<Family=F2,ID=INDI2,Father=0,Mother=0,Sex=2,Status=1>
##Sample=<Family=F3,ID=INDI3,Father=INDI1,Mother=INDI2,Sex=1,Status=1>
(...)

```





END_DOC
*/
import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
/**

BEGIN_DOC

This tools reads a pedigree file and inject it in the VCF header  


```
$ java -jar dist/vcfinjectpedigree.jar \
	-imih -imip -p input.ped \
	input.vcf.gz > out.vcf

$ grep Sample out.vcf
(...)
##Sample=<Family=F1,ID=INDI1,Father=0,Mother=0,Sex=1,Status=1>
##Sample=<Family=F2,ID=INDI2,Father=0,Mother=0,Sex=2,Status=1>
##Sample=<Family=F3,ID=INDI3,Father=INDI1,Mother=INDI2,Sex=1,Status=1>
(...)

```

END_DOC


 */
@Program(name="vcfinjectpedigree",
	description="Injects a pedigree (.ped) file in the VCF header",
	keywords={"vcf","pedigree","burden"})
public class VcfInjectPedigree
	extends Launcher
	{

private static final Logger LOG = Logger.build(VcfInjectPedigree.class).make();


@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
private File outputFile = null;


@Parameter(names={"-p","--pedigree"},description=Pedigree.OPT_DESCRIPTION)
private File pedigreeFile = null;

@Parameter(names={"-clean","--clean"},description="Remove all previous data about pedigree in the VCF header before adding the new one.")
private boolean cleanPreviousPedigree = false;

@Parameter(names={"-imih","--ignoreMissingInHeader"},description="Ignore errors if a sample is declared in the pedigree but is missing in the VCF header")
private boolean ignoreMissingInHeader = false;

@Parameter(names={"-imip","--ignoreMissingInPedigree"},description="Ignore errors if a sample is declared in the VCF header but is missing in the pedigree")
private boolean ignoreMissingInPedigree = false;

@Parameter(names={"-valid","--valid"},description="Ignore pedigree validation")
private boolean ignorePedigreeValidation = false;
	
	
	public VcfInjectPedigree()
		{
		}
	 
	
	
	@Override
	protected int doVcfToVcf(String inputName, VcfIterator in, VariantContextWriter out) {
		final VCFHeader header = in.getHeader();
		final List<String> sampleNames = Collections.unmodifiableList(header.getSampleNamesInOrder());

		
		try {
			LOG.info("reading "+this.pedigreeFile);
			final Pedigree pedigree = Pedigree.newParser().parse(this.pedigreeFile);
			
			if(!this.ignorePedigreeValidation)
				{
				pedigree.validate();
				}
			final Set<String> personsIds= new HashSet<>();
			for(final Pedigree.Person p:pedigree.getPersons()) {
				if(personsIds.contains(p.getId())) {
					LOG.error(
							"Pedigree contains two individual with ID :" +p.getId()
							+" This is not a pedigree error, but there is no mechanism to find them in the VCF header.");
					return -1;
				}
				personsIds.add(p.getId());
			}
			
			final Set<VCFHeaderLine> metaData = new HashSet<>(header.getMetaDataInInputOrder());
			
			/* remove previous pedigree */
			final Iterator<VCFHeaderLine> iter = metaData.iterator();
			while(iter.hasNext()) {
				final VCFHeaderLine hl = iter.next();
				if(Pedigree.VcfHeaderKey.equals( hl.getKey())) {
					if(this.cleanPreviousPedigree) {
					LOG.info("removing "+hl);
					iter.remove();
					}
				}
			}
			
			/* check already existing individual in vcf header */
			final Pedigree previousPedigree = Pedigree.newParser().parse(metaData);
			for(final Pedigree.Family fam:pedigree.getFamilies()) {
				final Pedigree.Family prevFam = previousPedigree.getFamilyById(fam.getId());
				if(prevFam==null) continue;
				for(final Pedigree.Person person:fam.getIndividuals()) {
					final Pedigree.Person prevPers = prevFam.getPersonById(person.getId());
					if(prevPers!=null) {
						LOG.error("VCFHeader from "+inputName+" already contains  sample "+person+
								". Use option -clean to removed it.");
						return -1;
					}
				}
			}
			
			/* check missing in header */
			for(final Pedigree.Person p: pedigree.getPersons()) {
				if(!sampleNames.contains(p.getId())) {
					if(!this.ignoreMissingInHeader) 
						{
						LOG.error("Sample "+p+" is declared in pedigree is missing in VCF header."+
								"Use option ----ignoreMissingInHeader  to disable this error.");
						return -1;
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
					if(!this.ignoreMissingInPedigree) 
					{
					LOG.error("VCF Sample "+sampleName+"  is missing in pedigree."+
							"Use option --ignoreMissingInPedigree to disable this error.");
					return -1;
					}
				else
					{
					LOG.warn("Sample "+sampleName+" is present in VCF header is missing in Pedigree.");
					}
				}
			}
			
			
			metaData.addAll(pedigree.toVCFHeaderLines());
			//super.addMetaData(metaData);
			
		
			final VCFHeader h2= new VCFHeader(metaData, sampleNames);
			final SAMSequenceDictionaryProgress progess=new SAMSequenceDictionaryProgress(header.getSequenceDictionary()).logger(LOG);
			
			out.writeHeader(h2);
			while(in.hasNext())
				{
				out.add(progess.watch(in.next()));
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
	public int doWork(List<String> args) {
		if(this.pedigreeFile==null || !this.pedigreeFile.exists()) {
			LOG.error("Undefined Pedigree file");
			return -1;
			}
		return doVcfToVcf(args, outputFile);
		}
	

	 	
	
	public static void main(String[] args)
		{
		new VcfInjectPedigree().instanceMainWithExit(args);
		}
	}
