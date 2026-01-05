/*
The MIT License (MIT)

Copyright (c) 2026 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.vcf2intervals;

import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bed.BedLine;
import com.github.lindenb.jvarkit.bed.BedLineReader;
import com.github.lindenb.jvarkit.bio.DistanceParser;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.dict.OrderChecker;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.variant.vcf.VCFReaderFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTextHeaderCodec;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CoordMath;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFReader;

/**
 
BEGIN_DOC

## Input

input is a VCF file or a VCF stream.
input must be sorted on chrom/pos.

## Example

```
$ java -jar dist/vcf2intervals.jar -N 5 src/test/resources/rotavirus_rf.vcf.gz  --min-distance 1
@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S1	SM:S1
@RG	ID:S2	SM:S2
@RG	ID:S3	SM:S3
@RG	ID:S4	SM:S4
@RG	ID:S5	SM:S5
@CO	vcf2intervals. compilation:20211112182935 githash:9b2ab03 htsjdk:2.24.1 date:20211112183125. cmd:-N 5 src/test/resources/rotavirus_rf.vcf.gz --min-distance 1
RF01	970	970	1	1
RF02	251	1965	5	1715
RF03	1221	2150	5	930
RF03	2201	2573	3	373
RF04	887	1860	5	974
RF04	1900	1920	2	21
RF05	41	1297	5	1257
RF05	1339	1339	1	1
RF06	517	1132	5	616
RF07	98	952	4	855
RF08	926	992	2	67
RF09	294	414	3	121
RF10	46	175	3	130
RF11	74	79	1	6
```

```
$ java -jar dist/vcf2intervals.jar --distance 300 --min-distance 0 src/test/resources/rotavirus_rf.vcf.gz  
@HD	VN:1.6	SO:coordinate
@SQ	SN:RF01	LN:3302
@SQ	SN:RF02	LN:2687
@SQ	SN:RF03	LN:2592
@SQ	SN:RF04	LN:2362
@SQ	SN:RF05	LN:1579
@SQ	SN:RF06	LN:1356
@SQ	SN:RF07	LN:1074
@SQ	SN:RF08	LN:1059
@SQ	SN:RF09	LN:1062
@SQ	SN:RF10	LN:751
@SQ	SN:RF11	LN:666
@RG	ID:S1	SM:S1
@RG	ID:S2	SM:S2
@RG	ID:S3	SM:S3
@RG	ID:S4	SM:S4
@RG	ID:S5	SM:S5
@CO	vcf2intervals. compilation:20211112182935 githash:9b2ab03 htsjdk:2.24.1 date:20211112183310. cmd:--distance 300 --min-distance 0 src/test/resources/rotavirus_rf.vcf.gz
RF01	970	970	1	1
RF02	251	251	1	1
RF02	578	877	2	300
RF02	1726	1965	2	240
RF03	1221	1242	2	22
RF03	1688	1708	2	21
RF03	2150	2315	3	166
RF03	2573	2573	1	1
RF04	887	991	2	105
RF04	1241	1262	2	22
RF04	1857	1920	3	64
RF05	41	41	1	1
RF05	499	795	2	297
RF05	879	879	1	1
RF05	1297	1339	2	43
RF06	517	695	4	179
RF06	1129	1132	1	4
RF07	98	225	2	128
RF07	684	952	2	269
RF08	926	992	2	67
RF09	294	414	3	121
RF10	46	175	3	130
RF11	74	79	1	6
```

END_DOC

*/

@Program(
	name="vcf2intervals",
	description="split a vcf to interval or bed for parallelization",
	keywords={"vcf","bed","interval"},
	creationDate="20211112",
	modificationDate="20221128",
	biostars= {9506628,9529137},
	jvarkit_amalgamion = true,
	menu="VCF Manipulation"
	)
public class VcfToIntervals extends Launcher
	{
	private static final Logger LOG=Logger.of(VcfToIntervals.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"--bed","--bed-output"},description="force BED format as output. (Default is '"+FileExtensions.INTERVAL_LIST+"')")
	private boolean force_bed_output = false;
	@Parameter(names={"-N","--variants","--n-variants"},description="number of variants per interval (or use option -D)")
	private long n_variants_per_interval = -1L;
	@Parameter(names={"-D","--distance"},description="min size of an interval (or use option -N). "+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int distance_per_interval = -1;
	@Parameter(names="--min-distance",description="extends the interval if the last variant is withing distance 'x' of the next interval. Ignore if negative."+DistanceParser.OPT_DESCRIPTION,converter=DistanceParser.StringConverter.class,splitter=NoSplitter.class)
	private int min_distance = -1;
	@Parameter(names={"--intervals","--bed-input"},description="Search for intervals for EACH record of the provided bed file. VCF path must be provided and indexed.")
	private Path bedIn = null;

	
	private void writeHeader(final PrintWriter pw,final SAMSequenceDictionary dict,final VCFHeader header ) {
		if(!this.force_bed_output) {
			final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
			final SAMFileHeader samHeader = new SAMFileHeader(dict);
			samHeader.setSortOrder(this.bedIn==null?SortOrder.coordinate:SortOrder.unknown);
			for(final String sn : header.getSampleNamesInOrder()) {
				final SAMReadGroupRecord rg = new SAMReadGroupRecord(sn);
				rg.setSample(sn);
				samHeader.addReadGroup(rg);
				}
			JVarkitVersion.getInstance().addMetaData(this, samHeader);
			codec.encode(pw, samHeader);
			}
		}
	
	private void apply(
			final PrintWriter pw,
			final PeekableIterator<VariantContext> iter,
			final OrderChecker<VariantContext> orderChecker,
			final Locatable optional_interval
			)
		{
		while(iter.hasNext()) {
			final VariantContext first = orderChecker.apply(iter.next());
			VariantContext last = first;
			long n_variants = 1;
			if(this.n_variants_per_interval>0L) {
				while(iter.hasNext() && n_variants < this.n_variants_per_interval) {
					if (!first.contigsMatch( iter.peek())) {
						break;
						}
					// consumme
					last = orderChecker.apply(iter.next());
					n_variants++;
					}
				}
			else {
				while(iter.hasNext()) {
					final VariantContext curr = iter.peek();
					if (!first.contigsMatch(curr)) {
						break;
						}
					if (CoordMath.getLength(first.getStart(), curr.getEnd()) > this.distance_per_interval) {
						break;
						}
					// consumme
					last = orderChecker.apply(iter.next());
					n_variants++;
					}
				}
			// next variant is just too close than the last one
			while(this.min_distance>=0 && iter.hasNext()) {
				final VariantContext curr = iter.peek();
				if(!last.withinDistanceOf(curr, this.min_distance)) break;
				// consumme
				last = orderChecker.apply(iter.next());
				n_variants++;
				}
			pw.print(first.getContig());
			pw.print("\t");
			pw.print(first.getStart()-(this.force_bed_output?1:0));
			pw.print("\t");
			pw.print(last.getEnd());
			pw.print("\t");
			pw.print(n_variants);
			pw.print("\t");
			pw.print(CoordMath.getLength(first.getStart(),last.getEnd()));
			pw.println();
			}// end while iter
		}
	
	@Override
	public int doWork(final List<String> args)
		{
		if(n_variants_per_interval>=0 && distance_per_interval>=0) {
			LOG.error("n-variants per interval and distance both defined");
			return -1;
			}
		else if(n_variants_per_interval<0 && distance_per_interval<0) {
			LOG.error("n-variants per interval or distance must be defined");
			return -1;
			}
		try
			{
			final String input = oneFileOrNull(args);
			if(StringUtils.isBlank(input) && bedIn!=null) {
				LOG.error("option --intervals doesn't work on stdin");
				return -1;
				}
			
			if(this.bedIn==null || StringUtils.isBlank(input)) {
				try(VCFIterator iter0  = super.openVCFIterator(input)) {
					final VCFHeader header =  iter0.getHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					final OrderChecker<VariantContext> orderChecker = new OrderChecker<>(dict, false);

					try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
						writeHeader(pw,dict,header);
						try(PeekableIterator<VariantContext> iter = new PeekableIterator<>(iter0)) {
							apply(pw,iter,orderChecker,null);
							}
						pw.flush();
						}//end writer
					}//end vcf reader
				}
			else
				{
				try(VCFReader reader = VCFReaderFactory.makeDefault().open(input,true)) {
					final VCFHeader header =  reader.getHeader();
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(header);
					final OrderChecker<VariantContext> orderChecker = new OrderChecker<>(dict, false);
					/* loop over each record in the bed */
					try(BedLineReader blr = new BedLineReader(this.bedIn)) {
						blr.setContigNameConverter(ContigNameConverter.fromOneDictionary(dict));
						try(PrintWriter pw = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
							writeHeader(pw,dict,header);
							while(blr.hasNext()) {
								final BedLine rec = blr.next();
								try(CloseableIterator<VariantContext> iter0 = reader.query(rec)) {
									try(PeekableIterator<VariantContext> iter = new PeekableIterator<>(iter0)) {
										apply(pw,iter,orderChecker,rec);
										}
									}
								}
							pw.flush();
							}//end pw
						}//end bed reader
					} // end vcf reader
				}
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}

	
	
	public static void main(final String[] args)
		{
		new VcfToIntervals().instanceMainWithExit(args);
		}

	}
