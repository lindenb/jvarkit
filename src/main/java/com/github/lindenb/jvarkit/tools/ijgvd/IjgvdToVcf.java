package com.github.lindenb.jvarkit.tools.ijgvd;
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



*/
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.AbstractIterator;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.samtools.ContigDictComparator;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.variant.vcf.VCFStandardHeaderLines;

/**

BEGIN_DOC

## DEPRECATED

Deprecated since data are now available as VCF

## DESCRIPTION

Integrative Japanese Genome Variation (iJGVD https://ijgvd.megabank.tohoku.ac.jp/ ) provides data of genomic variations obtained by whole-genome sequencing of Japanese individuals, who participate in the genome cohort study by ToMMo, IMM and other cohort projects in Japan.

> Rare variant discovery by deep whole-genome sequencing of 1,070 Japanese individuals, Nagasaki M, Yasuda J, Katsuoka F, Nariai N, Kojima K, Kawai Y, Yamaguchi-Kabat
a Y, Yokozawa J, Danjoh I, Saito S, Sato Y, Mimori T, Tsuda K, Saito R, Pan X, Nishikawa S, Ito S, Kuroki Y, Tanabe O, Fuse N, Kuriyama S, Kiyomoto H, Hozawa A, Minegi
shi N, Douglas Engel J, Kinoshita K, Kure S, Yaegashi N, ToMMo Japanese Reference Panel Project and Yamamoto M, Nat Commun, 21; 6:8018 (2015) 


##Example


```
$ ls ~/Downloads/chr* ~/Downloads/multiallelic_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip 
~/Downloads/chr10_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr11_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr1_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr12_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr13_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr14_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr15_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr16_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr17_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr18_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr19_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr20_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr21_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr2_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr22_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr3_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr4_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr5_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr6_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr7_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr8_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/chr9_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip
~/Downloads/multiallelic_20190717170045-bd410447a29ca7b9eb95ed92bf11b473.zip


$ java -jar dist/ijgv2vcf.jar -R ~/src/jvarkit-git/src/test/resources/human_b37.dict ~/Downloads/*.zip > out.vcf

$ java -jar dist/ijgv2vcf.jar -F -M -R ~/src/jvarkit-git/src/test/resources/human_b37.dict ~/Downloads/*.zip | bcftools view -O z -o ~/Downloads/3.5KJPN_tommo_2019071.vcf.gz
$ ls -lah ~/Downloads/3.5KJPN_tommo_2019071.vcf.gz
-rw-r--r-- 1 lindenb lindenb 514M juil. 17 15:04 ~/Downloads/3.5KJPN_tommo_2019071.vcf.gz
```

END_DOC

 */
@Program(name="ijgv2vcf",
	keywords={"vcf","jgvd","japan","tommo"},
	description="Convert zips of Integrative Japanese Genome Variation to VCF file.",
	creationDate="20190717",
	modificationDate="20190717",
	deprecatedMsg="Deprecated since data are now available as VCF")
public class IjgvdToVcf extends Launcher
	{
	private static final Logger LOG=Logger.build(IjgvdToVcf.class).make();
	
	private static final String REF_ALLELE_FREQ="REF_AF";
	private static final String REF_ALLELE_COUNT="REF_AC";
	private static final String ALT_ALLELE_FREQ="ALT_AF";
	private static final String ALT_ALLELE_COUNT="ALT_AC";
	private static final String TOTAL_ALLELES_COUNT="ALT_COUNT";
	private static final String N_SAMPLES="N_SAMPLES";
	private static final String FILTER1="FILTERED_IN_IJGVD";
	private static final String MULTIALLELIC="MULTIALLELIC_IN_IJGVD";
	
	

	
	private static final String IJGC_HEADERS[]=new String[] {
		"Chromosome",
		"Position_hg19",
		"rsSNP#",
		"Cytoband",
		"Variation_type",
		"Ref_allele",
		"Ref_allele_freq",
		"Ref_allele_count",
		"Alt_allele",
		"Alt_allele_freq",
		"Alt_allele_count",
		"Total_allele_count",
		"#Samples",
		"Gene_annotation",
		"Status",
		"Multi-allelic_info"
		};
    @Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
    private Path out =null;
    @Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
    private Path fai =null;
    @Parameter(names={"-F","--no-filtered"},description="ignore 'filtered' entries")
    private boolean skip_filtered= false;
    @Parameter(names={"-M","--no-multiallelic"},description="ignore 'multiallelic' entries")
    private boolean skip_multiallelic= false;
    @ParametersDelegate
    private WritingVariantsDelegate writingVariantsDelegate = new WritingVariantsDelegate();

    
    private ContigNameConverter ctgNameConverter = null;

    
    private class ZipIterator
    	extends AbstractIterator<VariantContext>
    	implements CloseableIterator<VariantContext>
    	{
    	final String entryName;
    	private BufferedReader br;
    	ZipIterator(final String entryName,InputStream in) throws IOException {
    		this.entryName = entryName;
    		this.br = new BufferedReader(new InputStreamReader(in, "UTF-8"));
    		}
    	@Override
    	protected VariantContext advance() {
    		try {
	    		for(;;)
	    			{
	    			final String line= this.br.readLine();
	    			if(line==null) return null;
	    			if(StringUtils.isBlank(line)) continue;
	    			if(line.startsWith("#")) continue;
	    			final String tokens[] = CharSplitter.TAB.split(line);
	    			if(line.startsWith("Chromosome")) {
	    				if(tokens.length>IJGC_HEADERS.length || IntStream.range(0, tokens.length).anyMatch(X->!tokens[X].equals(IJGC_HEADERS[X]))) {
	    					throw new RuntimeIOException("bad header in "+entryName+" got "+line+" but expected "+String.join("\t", IJGC_HEADERS));
	    					}
	    				continue;
	    				}
	    			final String ctg = ctgNameConverter.apply(tokens[0]);
	    			if(StringUtils.isBlank(ctg)) continue;
	    			final VariantContextBuilder vcb = new VariantContextBuilder();
	    			vcb.chr(ctg);
	    			
	    			final Allele ref= Allele.create(tokens[5],true);
	    			
	    			vcb.start(Integer.parseInt(tokens[1]));
	    			vcb.stop(Integer.parseInt(tokens[1])+tokens[5].length()-1);

	    			
	    			
	    			final List<Allele> alts= Arrays.stream(CharSplitter.COMMA.split(tokens[8])).map(A->Allele.create(A,false)).collect(Collectors.toList());
	    			
	    			final List<Allele> alleles = new ArrayList<>(alts.size()+1);
	    			alleles.add(ref);
	    			alleles.addAll(alts);
	    			vcb.alleles(alleles);
	    			
	    			if(!StringUtils.isBlank(tokens[2])) {
	    				vcb.id(tokens[2]);
	    				}
	    			
	    			vcb.attribute(REF_ALLELE_FREQ, Double.parseDouble(tokens[6]));
	    			vcb.attribute(REF_ALLELE_COUNT, Integer.parseInt(tokens[7]));
	    			vcb.attribute(ALT_ALLELE_FREQ, Double.parseDouble(tokens[9]));
	    			vcb.attribute(ALT_ALLELE_COUNT, Integer.parseInt(tokens[10]));
	    			vcb.attribute(TOTAL_ALLELES_COUNT, Integer.parseInt(tokens[11]));
	    			vcb.attribute(N_SAMPLES, Integer.parseInt(tokens[12]));
	    			
	    			if(alts.size()==1)
	    				{
	    				vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, Double.parseDouble(tokens[9]));
	    				vcb.attribute(VCFConstants.ALLELE_COUNT_KEY, Integer.parseInt(tokens[10]));
	    				}
	    			else
	    				{
	    				/** A:0.000281373|C:0.000140687 */
	    				final String tokens2[] = CharSplitter.PIPE.split(tokens[15]);
	    				vcb.attribute(VCFConstants.ALLELE_FREQUENCY_KEY, Arrays.stream(tokens2).mapToDouble(S->Double.parseDouble(S.substring(1+S.indexOf(":")))).toArray());
	    				}
	    			if(entryName.contains("filtered")) {
	    				vcb.filter(FILTER1);
	    				} 
	    			else if(alts.size()>1)
	    				{
	    				vcb.filter(MULTIALLELIC);
	    				}
	    			else
	    				{
	    				vcb.passFilters();
	    				}
	    			
	    			return vcb.make();
	    			}
    			}
    		catch(final IOException err)
    			{
    			throw new RuntimeIOException(err);
    			}
    		}
    	@Override
    	public void close() {
    		CloserUtil.close(this.br);
    		}
    	}

    @Override
    public int doWork(List<String> args) {
    	try {
    		final SAMSequenceDictionary dict= SequenceDictionaryUtils.extractRequired(this.fai);
    		this.ctgNameConverter = ContigNameConverter.fromOneDictionary(dict);
    		List<Path> zipPaths = IOUtils.unrollPaths(args);
			List<CloseableIterator<VariantContext>> iterators = new ArrayList<>(zipPaths.size()*2);
			for(final Path zipPath:zipPaths)
				{
				if(zipPath.getFileName().toString().endsWith(".tsv"))
					{
					String fname= zipPath.getFileName().toString();
					if( fname.endsWith("filtered.tsv") && skip_filtered)
						{
						continue;
						}
					if( fname.endsWith("_multiallelic.tsv") && skip_multiallelic)
						{
						continue;
						}
					final InputStream in =Files.newInputStream(zipPath);
					iterators.add(new ZipIterator(fname, in));
					continue;
					}
				for(int i=0;i<2;i++)
					{
					final ZipInputStream zin = new ZipInputStream(Files.newInputStream(zipPath));
					ZipEntry entry=null;
					ZipIterator zipIter = null;
					while((entry=zin.getNextEntry())!=null)
						{
						if(entry.getName().endsWith("filtered.tsv") && skip_filtered)
							{
							zin.closeEntry();
							continue;
							}
						if(entry.getName().endsWith("_multiallelic.tsv") && skip_multiallelic)
							{
							zin.closeEntry();
							continue;
							}
						if(i==0 && !(entry.getName().endsWith("passed.tsv")||entry.getName().endsWith("multiallelic.tsv"))) {
							zin.closeEntry();
							continue;
							}
						if(i==1 && !entry.getName().endsWith("filtered.tsv")) {
							zin.closeEntry();
							continue;
							}
						zipIter = new ZipIterator(entry.getName(),zin);
						break;
						}
					if(zipIter==null)
						{
						zin.close();
						}
					else
						{
						iterators.add(zipIter);
						}
					}
				}
			final ContigDictComparator contigDictComparator = new ContigDictComparator(dict);
			final Comparator<VariantContext> comparator = (A,B)->{
				int i = contigDictComparator.compare(A.getContig(), B.getContig());
				if(i!=0) return i;
				i= Integer.compare(A.getStart(), B.getStart());
				if(i!=0) return i;
				return A.getReference().compareTo(B.getReference());
			};
			final MergingIterator<VariantContext> iter = new MergingIterator<>(comparator, iterators);
			
			
			final VariantContextWriter vcw = writingVariantsDelegate.dictionary(dict).open(out);
			
			final Set<VCFHeaderLine> metaData = new HashSet<>();
			VCFStandardHeaderLines.addStandardInfoLines(metaData, true,
					VCFConstants.ALLELE_COUNT_KEY,
					VCFConstants.ALLELE_FREQUENCY_KEY,
					VCFConstants.ALLELE_NUMBER_KEY
					);
			
			metaData.add(new VCFInfoHeaderLine(REF_ALLELE_FREQ,1,VCFHeaderLineType.Float,"Ref Allele Freq."));
			metaData.add(new VCFInfoHeaderLine(ALT_ALLELE_FREQ,1,VCFHeaderLineType.Float,"Alt Allele Freq."));
			metaData.add(new VCFInfoHeaderLine(REF_ALLELE_COUNT,1,VCFHeaderLineType.Integer,"Ref Allele count"));
			metaData.add(new VCFInfoHeaderLine(ALT_ALLELE_COUNT,1,VCFHeaderLineType.Integer,"Alt Allele count."));
			metaData.add(new VCFInfoHeaderLine(TOTAL_ALLELES_COUNT,1,VCFHeaderLineType.Integer,"Total allele count"));
			metaData.add(new VCFInfoHeaderLine(N_SAMPLES,1,VCFHeaderLineType.Integer,"N samples"));
			metaData.add(new VCFFilterHeaderLine(MULTIALLELIC,"multiallelic"));
			metaData.add(new VCFFilterHeaderLine(FILTER1,"filtered in input"));
			
			
			VCFHeader header= new VCFHeader(metaData);
			header.setSequenceDictionary(dict);
			JVarkitVersion.getInstance().addMetaData(this, header);
			vcw.writeHeader(header);
			ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().dictionary(dict).logger(LOG).build();
			while(iter.hasNext()) {
				final VariantContext ctx = progress.apply(iter.next());
				vcw.add(ctx);
				}
			vcw.close();
			iter.close();
			progress.close();
			return 0;
		} catch (Exception e) {
			LOG.error(e);
			return -1;
			}
    	}
    
	public static void main(final String[] args) {
		new IjgvdToVcf().instanceMainWithExit(args);
		}
	}
