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


*/
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.picard.SAMSequenceDictionaryProgress;
import com.github.lindenb.jvarkit.util.vcf.TabixVcfFileReader;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import com.github.lindenb.jvarkit.util.vcf.VcfIterator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**

BEGIN_DOC


VCF files should be sorted using the same order as the sequence dictionary (see picard SortVcf).


### Example

list variants found with gatk AND samtools, keep the variants with http://www.sequenceontology.org/browser/current_release/term/SO:0001818 , remove variants found in a previous alignment (samtools or gatk)


```
#!/bin/bash

ls Samples | while read S
do
gunzip -c NEWALIGN/{S}.gatk.vcf.gz |\
        java -jar jvarkit-git/dist/vcffilterso.jar -A SO:0001818 |\
        java -jar jvarkit-git/dist/vcfin.jar NEWALIGN/{S}.samtools.vcf.gz |\
        java -jar jvarkit-git/dist/vcfin.jar -i OLDALIGN/{S}.samtools.vcf.gz |
        java -jar jvarkit-git/dist/vcfin.jar -i OLDALIGN/${S}.gatk.vcf.gz |
        grep -vE '^#' |
        awk -v S=${S} -F '      ' '{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n",S,$1,$2,$3,$4,$5,$8);}' 
done
```




#### Example 2

My list of bad variants is in the file 'bad.vcf'.



```
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	11167517	.	G	A	.	.	.
1	11167760	.	C	T	.	.	.
1	11168529	.	A	G	.	.	.
1	11169789	.	A	G	.	.	.
1	11174331	.	T	C	.	.	.
1	11174715	.	T	C	.	.	.
1	11180949	.	C	T	.	.	.
```



My main vcf file is 'input.vcf'.



```
##fileformat=VCFv4.2
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	11166480	.	C	T	.	.	.
1	11166541	.	G	C	.	.	.
1	11166577	.	C	T	.	.	.
1	11166713	.	T	C	.	.	.
1	11167146	.	G	A	.	.	.
1	11167158	.	C	T	.	.	.
1	11167270	.	G	T	.	.	.
1	11167517	.	G	T	.	.	.
1	11167627	.	G	C	.	.	.
1	11167760	.	C	T	.	.	.
1	11167829	.	C	T	.	.	.
1	11168529	.	A	G	.	.	.
1	11168769	.	CAAA	C	.	.	.
1	11169250	.	G	A	.	.	.
1	11169420	.	G	A	.	.	.
1	11169440	.	A	G	.	.	.
1	11169585	.	A	G	.	.	.
1	11169624	.	T	C	.	.	.
```



I want to put a FILTER in the variants if they are contained in VCF.



```
java -jar dist/vcfin.jar -A -fi InMyListOfBadVariants jeter2.vcf jeter1.vcf
```



output:



```
##fileformat=VCFv4.2
##FILTER=<ID=InMyListOfBadVariants,Description="Variant overlapping database.">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	11166480	.	C	T	.	.	.
1	11166541	.	G	C	.	.	.
1	11166577	.	C	T	.	.	.
1	11166713	.	T	C	.	.	.
1	11167146	.	G	A	.	.	.
1	11167158	.	C	T	.	.	.
1	11167270	.	G	T	.	.	.
1	11167517	.	G	T	.	.	.
1	11167627	.	G	C	.	.	.
1	11167760	.	C	T	.	InMyListOfBadVariants	.
1	11167829	.	C	T	.	.	.
1	11168529	.	A	G	.	InMyListOfBadVariants	.
1	11168769	.	CAAA	C	.	.	.
1	11169250	.	G	A	.	.	.
1	11169420	.	G	A	.	.	.
1	11169440	.	A	G	.	.	.
1	11169585	.	A	G	.	.	.
1	11169624	.	T	C	.	.	.

```


Please note that variant 1	11167517 is not flagged because is alternate allele is not contained in 'bad.vcf'




### History


 *  2015-02-24: rewritten. all files must be sorted: avoid to sort on disk. Support for tabix. Option -A
 *  2015-01-26: changed option '-v' to option '-i' (-v is for version)
 *  2014: Creation





END_DOC
*/

@Program(name="vcfin",
	description="Only prints variants that are contained/not contained into another VCF",
	keywords={"vcf","compare"},
	biostars={287815}
	)
public class VcfIn extends Launcher
	{
	private static final Logger LOG = Logger.build(VcfIn.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;
	@Parameter(names={"-i","--inverse"},description="Print variant that are not part of the VCF-database.")
	private boolean inverse = false;
	@Parameter(names={"-t","--tabix","--tribble","--indexed"},description="Database is indexed with tabix or tribble")
	private boolean databaseIsIndexed = false;
	@Parameter(names={"-A","--allalt"},description="ALL user ALT must be found in VCF-database ALT")
	private boolean userAltInDatabase = false;
	@Parameter(names={"-fi","--filterin"},description="Do not discard variant but add this FILTER if the variant is found in the database")
	private String filterIn = "";
	@Parameter(names={"-fo","--filterout"},description="Do not discard variant but add this FILTER if the variant is NOT found in the database")
	private String filterOut = "";

	public VcfIn()
		{
		}
		
	
	private boolean sameContext(			
			final VariantContext ctx1,
			final VariantContext ctx2
			)
		{
		return ctx1.getContig().equals(ctx2.getContig()) &&
				ctx1.getStart() == ctx2.getStart() &&
				ctx1.getEnd() == ctx2.getEnd() &&
				ctx1.getReference().equals(ctx2.getReference())
				;
		}
	
	private boolean allUserAltFoundInDatabase(
			final VariantContext userVariants,
			final VariantContext databaseVariants
			)
		{
		if(!this.userAltInDatabase) return true;
		final Set<Allele> user_alts=new HashSet<Allele>(userVariants.getAlternateAlleles());
		user_alts.removeAll(databaseVariants.getAlternateAlleles());
		return user_alts.isEmpty();
		}
	
	@Override
	protected VCFHeader addMetaData(final VCFHeader header) {
		if(!this.filterIn.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(this.filterIn,
					"Variant overlapping database."));
			}
		if(!this.filterOut.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(this.filterOut,
					"Variant non overlapping database."));
			}
		return super.addMetaData(header);
		}
	
	private void addVariant(final VariantContextWriter w,final VariantContext ctx,boolean keep)
		{
		if(this.inverse) keep=!keep;
		if(!this.filterIn.isEmpty())
			{
			if(keep){
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(this.filterIn);
				w.add(vcb.make());
				}
			else
				{
				w.add(ctx);
				}
			}
		else  if(!this.filterOut.isEmpty()) {
			if(keep){
				w.add(ctx);
				}
			else
				{
				final VariantContextBuilder vcb=new VariantContextBuilder(ctx);
				vcb.filter(this.filterOut);
				w.add(vcb.make());
				}
			}
		else
			{
			if(keep) {
				w.add(ctx);
				} else
				{
					/* don't print */
				}
			}
		}
	
	private int scanFileSorted(
			final VariantContextWriter vcw,
			final String databaseVcfUri,
			final VcfIterator userVcfIn
			)
		{
		EqualRangeVcfIterator equalRangeDbIter=null;
		EqualRangeIterator<VariantContext> equalRangeUserVcf = null;
		try
			{
			final VCFHeader header = new VCFHeader(userVcfIn.getHeader());
			final SAMSequenceDictionary userVcfDict = header.getSequenceDictionary();
			/// NO need if(dict1==null)
			if(userVcfDict==null)
				{
				LOG.error(JvarkitException.VcfDictionaryMissing.getMessage("user file"));
				return -1;
				}
			final Comparator<VariantContext> userVcfComparator =
					VCFUtils.createTidPosComparator(userVcfDict)
					;
			equalRangeDbIter = new EqualRangeVcfIterator(
					VCFUtils.createVcfIterator(databaseVcfUri),userVcfComparator);

			this.addMetaData(header);
			vcw.writeHeader(header);
			final SAMSequenceDictionaryProgress progress = new SAMSequenceDictionaryProgress(userVcfDict).logger(LOG);
			
			equalRangeUserVcf = new EqualRangeIterator<>(userVcfIn, userVcfComparator);
			
			while(equalRangeUserVcf.hasNext())
				{
				final List<VariantContext> ctxList = equalRangeUserVcf.next();
				progress.watch(ctxList.get(0));
				
				//fill both contextes
				final List<VariantContext> dbContexes = new ArrayList<VariantContext>(equalRangeDbIter.next(ctxList.get(0)));
				
				for(final VariantContext userCtx:ctxList)
					{
					boolean keep = dbContexes.stream().
							filter(V->sameContext(userCtx,V)).
							anyMatch(V->allUserAltFoundInDatabase(userCtx,V))
							;
					addVariant(vcw,userCtx,keep);
					}
				if(vcw.checkError()) break;
				}
			equalRangeUserVcf.close();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(equalRangeDbIter);
			CloserUtil.close(userVcfIn);
			CloserUtil.close(vcw);
			}
		}

	private int scanUsingTabix(
			final VariantContextWriter vcw,
			final String databaseFile,
			final VcfIterator in2
			)
		{
		VCFFileReader tabix=null;
		try
			{
			tabix =  new VCFFileReader(new File(databaseFile),true);
			final VCFHeader header1= new VCFHeader(in2.getHeader());
			this.addMetaData(header1);
			vcw.writeHeader(header1);
			
			final SAMSequenceDictionaryProgress progress=new SAMSequenceDictionaryProgress(header1.getSequenceDictionary()).logger(LOG);
			
			while(in2.hasNext() && !vcw.checkError())
				{
				final VariantContext userCtx= progress.watch(in2.next());
				final CloseableIterator<VariantContext> iter= tabix.query(
						userCtx.getContig(),
						Math.max(1,userCtx.getStart()-1),
						userCtx.getEnd()+1
						);
				boolean keep=false;
				while(iter.hasNext())
					{
					final VariantContext dbctx= iter.next();
					if(!sameContext(userCtx,dbctx)) continue;
					if(!allUserAltFoundInDatabase(userCtx, dbctx)) continue;
					keep=true;
					break;
					}
				iter.close();
				
				addVariant(vcw,userCtx,keep);
				if(vcw.checkError()) break;
				}
			progress.finish();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tabix);
			CloserUtil.close(in2);
			}
		}
	@Override
	public int doWork(final List<String> args) {
		if(!this.filterIn.isEmpty() && !this.filterOut.isEmpty()) {
			 LOG.error("Option filterIn/filterOut both defined.");
			 return -1;
		}
		if(this.inverse && (!this.filterIn.isEmpty() || !this.filterOut.isEmpty())) {
			 LOG.error("Option inverse cannot be used when Option filterin/filterou is defined.");
			 return -1;
		}
		
		String databaseVcfUri;
		String userVcfUri;
		if(args.size()==1)
			{
			databaseVcfUri = args.get(0);
			userVcfUri =null;
			}
		else if(args.size()==2)
			{
			databaseVcfUri = args.get(0);
			userVcfUri = args.get(1);
			}
		else
			{
			LOG.error("illegal number of arguments");
			return -1;
			}

		VariantContextWriter w=null;
		VcfIterator in=null;
		try {
			in = (userVcfUri==null?
					VCFUtils.createVcfIteratorFromInputStream(stdin()):
					VCFUtils.createVcfIterator(userVcfUri)
					);
			w= super.openVariantContextWriter(outputFile);
			if(this.databaseIsIndexed)
				{
				return this.scanUsingTabix(w,databaseVcfUri, in);
				}
			else
				{
				return this.scanFileSorted(w,databaseVcfUri, in);
				}
			} catch (final Exception err) {
				LOG.error(err);
				return -1;
			} finally
			{
			CloserUtil.close(in);
			CloserUtil.close(w);
			}
		}

	public static void main(final String[] args) {
		new VcfIn().instanceMainWithExit(args);
	}
	}
