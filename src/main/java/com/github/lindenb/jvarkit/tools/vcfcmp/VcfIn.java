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


*/
package com.github.lindenb.jvarkit.tools.vcfcmp;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.iterator.EqualRangeIterator;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;
import htsjdk.variant.vcf.VCFIterator;

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
import htsjdk.samtools.util.StringUtil;

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
        java -jar jvarkit-git/dist/vcfin.jar -D NEWALIGN/{S}.samtools.vcf.gz |\
        java -jar jvarkit-git/dist/vcfin.jar -i -D OLDALIGN/{S}.samtools.vcf.gz |
        java -jar jvarkit-git/dist/vcfin.jar -i -D OLDALIGN/${S}.gatk.vcf.gz |
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
java -jar dist/vcfin.jar -A -fi InMyListOfBadVariants -D jeter2.vcf jeter1.vcf
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

 *  2018-08-31: database must be specified with '-D'. ContigName conversion applied for tabix/tribble data.
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
	@Parameter(names={"-t","--tabix","--tribble","--indexed"},description=
			"Database is indexed with tabix or tribble. I will use random access to get the data. Otherwise, VCF are assumed to be sorted the same way."
			+ " [20180731] I will try to convert the contig names e.g: 'chr1' -> '1'")
	private boolean databaseIsIndexed = false;
	@Parameter(names={"-A","--allalt"},description="ALL user ALT must be found in VCF-database ALT")
	private boolean userAltInDatabase = false;
	@Parameter(names={"-fi","--filterin"},description="Do not discard variant but add this FILTER if the variant is found in the database(s)")
	private String filterIn = "";
	@Parameter(names={"-fo","--filterout"},description="Do not discard variant but add this FILTER if the variant is NOT found in the database(s)")
	private String filterOut = "";
	@Parameter(names={"-D","--database"},description="external database uri",required=true)
	private List<String> externalDatabaseList = new ArrayList<>();
	
	// http://gnomad.broadinstitute.org/dbsnp/rs11361742
	// gnomad www: 19:45454285 TA / T 
	// gnomad www: 19:45454285 TAA/ T
	// my data: 19	45454285	TA	T
	// gnomad-vcf: 19	45454285	TAA	TA,T,TAAA

	@Parameter(names={"-cp","--only-contig-pos"},description=
			"Two variants are the same if they have the same CONTIG/POS. Do NOT Look at REF or ALTS. "
			+ "Motivation: e.g  http://gnomad.broadinstitute.org/dbsnp/rs11361742 "
			+ "in gnomad VCF: [19	45454285 TAA	TA,T,TAAA ] "
			+ "in my VCF: [19	45454285	TA	T]. "
			+ "Only works with --tabix option")
	private boolean only_contig_pos = false;
	@Parameter(names="-m",description="min number of equivalent variants found in database(s), inclusive. With --tabix mode only. .eg: '2': the user variant must be found in at least 2 VCF database.")
	private int minCountInclusive=1;
	@Parameter(names="-M",description=" max number of equivalent variants found in database(s). -1 : no limit. sinclusive.With --tabix mode only. .eg: '3': the user variant must be found in 3 or less VCF database.")
	private int maxCountInclusive= -1;


	public VcfIn()
		{
		}
		
	/** if dict are different , ignore contig*/
	private boolean sameContextIgnoreContig(			
			final VariantContext ctx1,
			final VariantContext ctx2
			)
		{
		return ctx1.getStart() == ctx2.getStart() &&
				ctx1.getEnd() == ctx2.getEnd() &&
				ctx1.getReference().equals(ctx2.getReference())
				;
		}

	
	private boolean sameContext(			
			final VariantContext ctx1,
			final VariantContext ctx2
			)
		{
		return ctx1.getContig().equals(ctx2.getContig()) &&
				sameContextIgnoreContig(ctx1,ctx2)
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
					"Variant overlapping database(s). " + 
			String.join(" ",this.externalDatabaseList)));
			}
		if(!this.filterOut.isEmpty()) {
			header.addMetaDataLine(new VCFFilterHeaderLine(this.filterOut,
					"Variant non overlapping database(s)." +
			String.join(" ",this.externalDatabaseList)));
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
				if(ctx.isFiltered()) {
					w.add(ctx);
					}
				else
					{
					w.add(new VariantContextBuilder(ctx).passFilters().make());
					}
				}
			}
		else  if(!this.filterOut.isEmpty()) {
			if(keep){
				if(ctx.isFiltered()) {
					w.add(ctx);
					}
				else
					{
					w.add(new VariantContextBuilder(ctx).passFilters().make());
					}
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
			final VCFIterator userVcfIn,
			String externalDatabaseURI
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
					VCFUtils.createVCFIterator(externalDatabaseURI),userVcfComparator);

			this.addMetaData(header);
			vcw.writeHeader(header);
			
			final ProgressFactory.Watcher<VariantContext> progress= 
					ProgressFactory.newInstance().
						dictionary(header).
						logger(LOG).
						build();
			
			equalRangeUserVcf = new EqualRangeIterator<>(userVcfIn, userVcfComparator);
			
			while(equalRangeUserVcf.hasNext())
				{
				final List<VariantContext> ctxList = equalRangeUserVcf.next();
				progress.apply(ctxList.get(0));
				
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
			progress.close();
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
	
	private class TabixVcf
		implements Closeable
		{
		final File file;
		final VCFFileReader vcfFileReader;
		final VCFHeader header;
		final SAMSequenceDictionary dict;
		final ContigNameConverter contigNameConverter;
		TabixVcf(SAMSequenceDictionary dictIn,final String file) {
			this.file = new File(file);
			this.vcfFileReader = new VCFFileReader(this.file,true);
			this.header = this.vcfFileReader.getFileHeader();
			this.dict = this.header.getSequenceDictionary();
			
			if(dictIn!=null && this.dict!=null)
				{
				contigNameConverter = ContigNameConverter.fromDictionaries(dictIn, this.dict);
				}
			else if(this.dict!=null)
				{
				contigNameConverter = ContigNameConverter.fromOneDictionary(this.dict);
				}
			else
				{
				contigNameConverter = ContigNameConverter.getIdentity();
				}
			}
		
		@Override
		public void close() throws IOException {
			CloserUtil.close(vcfFileReader);
			}
		}

	private int scanUsingTabix(
			final VariantContextWriter vcw,
			final VCFIterator vcfIn,
			final List<String> databasesPaths
			)
		{
		final List<TabixVcf> tabixList = new ArrayList<>(databasesPaths.size());
		try
			{
			final VCFHeader header1= new VCFHeader(vcfIn.getHeader());
			final SAMSequenceDictionary dictIn = header1.getSequenceDictionary();
			

			
			for(final String path: databasesPaths)
				{
				tabixList.add(new TabixVcf(dictIn,path));
				}
			
			final VCFHeader header2 = new VCFHeader(header1);
			this.addMetaData(header2);
			vcw.writeHeader(header2);
			
			final ProgressFactory.Watcher<VariantContext> progress= 
					ProgressFactory.newInstance().
						dictionary(dictIn).
						logger(LOG).
						build();
			
			while(vcfIn.hasNext())
				{
				final VariantContext userCtx= progress.apply(vcfIn.next());
				
				int number_of_time_ctx_was_found = 0;
				for(final TabixVcf tabix:tabixList)  {
					final String newContigName = tabix.contigNameConverter.apply(userCtx.getContig());
					if(StringUtil.isBlank(newContigName))  continue;
					final CloseableIterator<VariantContext> iter= tabix.vcfFileReader.query(
							newContigName,
							Math.max(1,userCtx.getStart()-1),
							userCtx.getEnd()+1
							);
					final Predicate<VariantContext> acceptVariant = (V)->{
						if(!V.getContig().equals(newContigName)) return false;
						if(only_contig_pos && V.getStart()==userCtx.getStart()) return true;
						if(!sameContextIgnoreContig(userCtx,V)) return false;
						if(!allUserAltFoundInDatabase(userCtx, V)) return false;
						return true;
						};
					
					while(iter.hasNext())
						{
						final VariantContext dbctx= iter.next();
						if(!acceptVariant.test(dbctx)) continue;
						number_of_time_ctx_was_found++;
						break;
						}
					iter.close();
					if( number_of_time_ctx_was_found >= minCountInclusive) break;
					}
				final boolean keep = number_of_time_ctx_was_found >= this.minCountInclusive &&
						(this.maxCountInclusive<0 || number_of_time_ctx_was_found <= this.maxCountInclusive)
						;
				addVariant(vcw,userCtx,keep);
				if(vcw.checkError()) break;
				}
			progress.close();
			return RETURN_OK;
			}
		catch(final Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(tabixList);
			CloserUtil.close(vcfIn);
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
		
		if(args.size()==2) {
			LOG.error("cmd-line syntax has changed. Database must be specified with --database/-D");
			return -1;
			}
		
		
		
		
		
		
		final String userVcfUri = oneFileOrNull(args);
		VariantContextWriter w=null;
		VCFIterator in=null;
		try {
			in = (userVcfUri==null?
					VCFUtils.createVCFIteratorFromInputStream(stdin()):
					VCFUtils.createVCFIterator(userVcfUri)
					);
			
			final List<String> externalDatabaseUris = 
					this.externalDatabaseList.size() ==1 && this.externalDatabaseList.get(0).endsWith(".list") ?
					Files.lines(Paths.get(this.externalDatabaseList.get(0))).filter(L->!StringUtil.isBlank(L)).collect(Collectors.toList())
					:
					this.externalDatabaseList
					;
			
			if(externalDatabaseUris.isEmpty()) {
				LOG.error("Empty external database URI");
				return -1;
				}
			
			w= super.openVariantContextWriter(outputFile);
			if(this.databaseIsIndexed)
				{
				return this.scanUsingTabix(w,in,externalDatabaseUris);
				}
			else
				{
				if(externalDatabaseUris.size()!=1) {
					LOG.info("When files are sorted, I can only use one URI.");
					return -1;
					}
				if(this.only_contig_pos)
					{
					LOG.error("--only-contig-pos only  works with --tabix flag ");
					return -1;
					}
				return this.scanFileSorted(w,in,externalDatabaseUris.get(0));
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
