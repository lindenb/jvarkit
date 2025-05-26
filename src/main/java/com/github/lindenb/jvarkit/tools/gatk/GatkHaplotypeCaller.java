/*
The MIT License (MIT)
Copyright (c) 2025 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gatk;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.Stream;


import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.gatk.Gatk4Proxy;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;
/**
BEGIN_DOC

This tool is a wrapper over gatk4 HaplotypeCaller

Not All versions of jvarkit will work for this tool. Because it needs to
be compiled against the java code from gatk using a command like:

```
./gradlew jvarkit -Dgatk4.local.jar=/path/to/gatk/gatk-package-4.*-local.jar
```

Furthermore, don't trust the automatic generated documentation,  the tool must be invoked using the following syntax:

```
java -cp /path/to/jvarkit.jar:/path/to/gatk-package-4.xxx-local.jar \
	com.github.lindenb.jvarkit.tools.gatk.GatkHaplotypeCaller \
	-R src/test/resources/rotavirus_rf.fa \
	--bed input.bed \
	src/test/resources/S*.bam
```

This tool:

* call each BAM into a .vcf.gz using haplotype caller
  * convert the dictionary/chromosomes if needed
  * convert the sample name if the same sample if present more than once
* group g.vcf.gz files by the sqrt(number-of-gvcf-files)
* combine each group of gvcf files
* genotypegvcfs for the final gvcf file

END_DOC

**/
@Program(name="gatkhc",
	description="Wrapper for GATK HaplotypeCaller",
	keywords={"gatk","vcf","bam"},
	modificationDate="20240625",
	creationDate="20240625",
	generate_doc = true,
	jvarkit_amalgamion = true
	)
public class GatkHaplotypeCaller extends AbstractGatkTool {
	private static final Logger LOG = Logger.of(GatkHaplotypeCaller.class);

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path referencePath = null;
	@Parameter(names={"--references"},description="Other references. If a reference is different from the main reference, the contigs of a  GVCF file will be converted (e.g: 1 -> chr1) to the main reference dictionary. ")
	private List<String> references = new ArrayList<>();
	@Parameter(names={"-L","-bed","--bed"},description="restrict to bed",required=true)
	private Path bedPath = null;
	@Parameter(names={"-dbsnp","--dbsnp"},description="path to dbsnp")
	private Path dbsnp = null;
	@Parameter(names={"--tmp","--tmp-dir"},description="temporary directory")
	private Path _baseTmpDir = IOUtils.getDefaultTempDir();
	@Parameter(names={"--mapq"},description="mapping quality")
	private int mapq=10;
	
	private Path tmpDir=null;

	private static class Reference {
		final Path path;
		final SAMSequenceDictionary dict;
		Reference(Path path) {
			this.path= path;
			this.dict=  SequenceDictionaryUtils.extractRequired(path);
			}
		final UnaryOperator<String> getContigConverter() {
			return ContigNameConverter.fromOneDictionary(dict);
			}
		@Override
		public boolean equals(Object obj) {
			if(obj==this) return true;
			if(obj==null || !(obj instanceof Reference)) return false;
			return areSequenceDictionariesEqual(Reference.class.cast(obj));
			}
		@Override
		public int hashCode() {
			return dict.hashCode();
			}
		
		boolean areSequenceDictionariesEqual(final Reference o) {
			return SequenceUtil.areSequenceDictionariesEqual(this.dict, o.dict);
			}
		
		@Override
		public String toString() {
			return path.toString();
			}
		}
	
	private Reference reference = null;
	
	private static class SamInput {
		String path;
		int index=-1;
		String srcSample;
		String newSample;
		Reference reference;
		@Override
		public String toString() {
			return path;
			}
		
		}
	

	
	

	private void fill_cmd(List<String> cmd,final Path bed,final Reference reference) {
		cmd.add("-R");
		cmd.add(reference.path.toString());
		//cmd.add("--verbosity");
		//cmd.add("ERROR");
		if(bed!=null) {
			cmd.add("-L");
			cmd.add(bed.toString());
			}
		cmd.add("--tmp-dir");
		cmd.add(this.tmpDir.toString());
		
		}

	private Path invoke_gvcf(final SamInput sampInput,final Path bed) throws Exception {
		final String basename=String.format("hc.%03d.g" + FileExtensions.COMPRESSED_VCF, sampInput.index);
		final Path gvcf =  this.tmpDir.resolve(basename);
		final Path bed0;
		
		if(sampInput.reference.equals(this.reference)) {
			bed0 = bed;
			}
		else
			{
			LOG.info("converting bed file "+bed);
			bed0  = Files.createTempFile(this.tmpDir, "tmp", FileExtensions.BED);
			final UnaryOperator<String> contigConverter = sampInput.reference.getContigConverter();
			boolean ok=false;
			try(BufferedReader br=IOUtils.openPathForBufferedReading(bed)) {
				try(BufferedWriter pw = Files.newBufferedWriter(bed0)) {
					String line;
					while((line=br.readLine())!=null) {
						String[] tokens= CharSplitter.TAB.split(line);
						final String ctg = contigConverter.apply(tokens[0]);
						if(StringUtils.isBlank(ctg)) continue;
						tokens[0]=ctg;
						pw.append( String.join("\t", tokens));
						pw.append("\n");
						ok=true;
						}
					pw.flush();
					}
				}
			if(!ok) {
				final String smallContig = sampInput.reference.dict.getSequences().
						stream().
						sorted((A,B)->Integer.compare(B.getLengthOnReference(), A.getLengthOnReference())).
						map(it->it.getContig()).
						findFirst().
						get();
				LOG.info("no compatible contif found for new reference in "+bed+". creating a ~empty bed with contig "+ smallContig);

				try(BufferedWriter pw = Files.newBufferedWriter(bed0)) {
					pw.append( smallContig+"\t0\t1\n");
					pw.flush();
					}
				}
			}
		
		final List<String> cmd= new ArrayList<>();
		cmd.add("HaplotypeCaller");
		fill_cmd(cmd,bed0,sampInput.reference);
		cmd.add("-I");
		cmd.add(sampInput.path);
		cmd.add("-O");
		cmd.add(gvcf.toString());
		cmd.add("-ERC");
		cmd.add("GVCF");
		if(this.mapq>0) {
			cmd.add("--minimum-mapping-quality");
			cmd.add(String.valueOf(this.mapq));
			}
		if(sampInput.reference.equals(this.reference) && this.dbsnp!=null) {
			cmd.add("--dbsnp");
			cmd.add(String.valueOf(this.dbsnp));
			}
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-G");cmd.add("StandardHCAnnotation");
		
		execute(cmd);
		
		/* delete converted bed file if needed */
		if(!sampInput.reference.equals(this.reference)) {
			FilesDelete(bed0);
			}
		
		
		if(!sampInput.reference.equals(this.reference)) {
			final UnaryOperator<String> contigConverter = this.reference.getContigConverter();
			
			
			final Path rename= Files.createTempFile("tmp", FileExtensions.COMPRESSED_VCF);
			try(VCFFileReader r=new VCFFileReader(gvcf,true)) {
				final VCFHeader header= r.getFileHeader();
				
				final VCFHeader header2 = new VCFHeader(header);
				header2.setSequenceDictionary(this.reference.dict);
				
				final  VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder().
						setOutputPath(rename).
						setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF).
						unsetOption(Options.INDEX_ON_THE_FLY).
						setCreateMD5(false).
						setReferenceDictionary(this.reference.dict).
						clearIndexCreator();
				try(VariantContextWriter vcw=vcwb.build()) {
					vcw.writeHeader(header2);
					try(CloseableIterator<VariantContext> iter=r.iterator()) {
						while(iter.hasNext()) {
							final VariantContext vc = iter.next();
							final String ctg = contigConverter.apply(vc.getContig());
							if(StringUtils.isBlank(ctg)) continue;
							if(ctg.equals(vc.getContig())) {
								vcw.add(vc);
								}
							else
								{
								vcw.add(new VariantContextBuilder(vc).chr(ctg).make());
								}
							}
						}
					}
				}
			
			FilesDelete(gvcf);
			FilesDelete(super.indexFor(gvcf));
			
			final Path rename2 = Files.createTempFile("tmp", FileExtensions.COMPRESSED_VCF);
			cmd.clear();
			cmd.add("SortVcf");
			cmd.add("--INPUT");
			cmd.add(rename.toString());
			cmd.add("--OUTPUT");
			cmd.add(rename2.toString());			
			cmd.add("--TMP_DIR");
			cmd.add(this.tmpDir.toString());			
			execute(cmd);
			
			
			FilesDelete(rename);
			indexFeatureFile(rename2);
			

			Files.move(rename2, gvcf);
			Files.move(indexFor(rename2),indexFor(gvcf));
			}
		
		if(!sampInput.srcSample.equals(sampInput.newSample)) {
			final Path rename= Files.createTempFile("tmp", FileExtensions.COMPRESSED_VCF);
			cmd.clear();
			cmd.add("RenameSampleInVcf");
			cmd.add("--INPUT");
			cmd.add(gvcf.toString());
			cmd.add("--OUTPUT");
			cmd.add(rename.toString());
			cmd.add("-NEW_SAMPLE_NAME");
			cmd.add(sampInput.newSample);
			cmd.add("--CREATE_INDEX");
			cmd.add("false");
			cmd.add("--TMP_DIR");
			cmd.add(this.tmpDir.toString());	
			cmd.add("--VALIDATION_STRINGENCY");
			cmd.add("LENIENT");
			
			execute(cmd);
			
			
			indexFeatureFile(rename);

			
			FilesDelete(gvcf);
			FilesDelete(super.indexFor(gvcf));
			Files.move(rename, gvcf);
			Files.move(indexFor(rename),indexFor(gvcf));
			}
		
		return gvcf;
		}

	private Path invoke_combine(int idx,final List<Path> gvcfs,final Path bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("combine.%03d.g" + FileExtensions.COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("CombineGVCFs");
		fill_cmd(cmd,bed,this.reference);
		for(Path gvcf:gvcfs) {
			cmd.add("-V");
			cmd.add(gvcf.toString());
			}
		cmd.add("-O");
		cmd.add(outPath.toString());
		if(this.dbsnp!=null) {
			cmd.add("--dbsnp");
			cmd.add(String.valueOf(this.dbsnp));
			}
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");

		
		execute(cmd);
		return outPath;
		}
	
	private void invoke_genotype(final Path gvcf,final Path bed,final Path outPath) throws Exception {
		final List<String> cmd= new ArrayList<>();
		final Path vcfOut;
		if(outPath==null) {
			vcfOut = Files.createTempFile(this.tmpDir,"hc.",".vcf.gz");
			}
		else
			{
			vcfOut=outPath;
			}
		cmd.add("GenotypeGVCFs");
		fill_cmd(cmd,bed,this.reference);
		cmd.add("-V");
		cmd.add(gvcf.toString());
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-O");
		cmd.add(vcfOut.toString());
		if(this.dbsnp!=null) {
			cmd.add("--dbsnp");
			cmd.add(String.valueOf(this.dbsnp));
			}
		execute(cmd);
		if(outPath==null) {
			try(VCFIterator r=new VCFIteratorBuilder().open(vcfOut)) {
				try(VariantContextWriter w= VCFUtils.createVariantContextWriterToOutputStream(stdout())) {
					w.writeHeader(r.getHeader());
					VCFUtils.copyVariantsTo(r, w);
					}
				}
			FilesDelete(vcfOut);
			FilesDeleteIfExists(indexFor(vcfOut));
			}
		}

	private void indexFeatureFile(final Path f) throws Exception {
		/** index VCF */
		final List<String> cmd=new ArrayList<>();
		cmd.add("IndexFeatureFile");
		cmd.add("--input");
		cmd.add(f.toString());
		execute(cmd);
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
		}
	
	@Override
	public int doWork(List<String> args) {
	try {
		if(getGatkEngine()==null) {
			System.err.println(Gatk4Proxy.MSG_COMPILATION);
			return -1;
			}
		
		IOUtil.assertDirectoryIsWritable(this._baseTmpDir);
		
		this.reference = new Reference(this.referencePath);
		
		final List<Reference> allRefs = 
				Stream.concat(
						Collections.singletonList(this.reference).stream(),
						IOUtils.unrollPaths(this.references).stream().map(P->new Reference(P))
						).collect(Collectors.toList());
					
		
		LOG.info("scan bams...");
		final List<SamInput> bams= new ArrayList<>();
		final Set<String> samples= new HashSet<>();
		for(String p: IOUtils.unrollStrings(args)) {
			final SamInput si = new SamInput();
			si.path = p;
			si.index = bams.size();
			final SamReaderFactory srf=SamReaderFactory.make();
			srf.validationStringency(ValidationStringency.LENIENT);
			final SAMFileHeader hdr;

			try(SamReader sfr=srf.open(SamInputResource.of(p))) {
				hdr=sfr.getFileHeader();
				}
			
			final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(hdr);
			
			si.reference = allRefs.stream().filter(R->SequenceUtil.areSequenceDictionariesEqual(dict, R.dict)).findFirst().orElseThrow(()->new IllegalArgumentException("Cannot find reference dictionary for BAM: p"));
			si.srcSample  = hdr.getReadGroups().stream().
					map(RG->RG.getSample()).
					filter(SN->!StringUtils.isBlank(SN)).
					findFirst().
					orElse(	
						IOUtils.getFilenameWithoutCommonSuffixes(Paths.get(p))
						);
			
			si.newSample = si.srcSample;
			int n=0;
			while(samples.contains(si.newSample)) {
				si.newSample = si.srcSample+":"+(++n);
				}
			if(!si.newSample.equals(si.srcSample)) {
				LOG.warn("Sample "+si.srcSample+" will be renamed "+si.newSample+" for "+p);
				}
			
			samples.add(si.srcSample);
			samples.add(si.newSample);
			bams.add(si);
			}
		
		if (bams.isEmpty()) {
        	LOG.error("No Bam was provided");
        	return -1;
        	}
		
	
		
		this.tmpDir = Files.createTempDirectory(this._baseTmpDir,"tmp");
		LOG.info("tmp directory will be "+this.tmpDir);
		
		final List<Path> gvcfs_list = new ArrayList<>(bams.size());
		for(int i=0;i< bams.size();i++) {
			LOG.info("("+(i+1)+"/"+bams.size()+" "+bams.get(i));
			final Path g_vcf_gz = invoke_gvcf(bams.get(i),this.bedPath);
			gvcfs_list.add(g_vcf_gz);
			}
			
		final int sqrt = Math.max(20, (int)Math.sqrt(gvcfs_list.size()));
		LOG.info("gvcfs will be pooled into clusters of "+sqrt+" GVCF files");

		final List<Path> combined0_list = new ArrayList<>(sqrt);
		int i=0;
		while(!gvcfs_list.isEmpty()) {
			final List<Path> L = new ArrayList<>(sqrt);
			while(!gvcfs_list.isEmpty() && L.size() < sqrt) {
				L.add(gvcfs_list.remove(0));
				}
			if(L.size()==1) {
				combined0_list.add(L.get(0));
				}
			else {
				final Path combined = invoke_combine(i,L, this.bedPath);
				combined0_list.add(combined);
				for(Path p:L) {
					FilesDelete(p);
					FilesDelete(indexFor(p));
					}
				}
			i++;
			}
		Path to_genotype;
		if(combined0_list.size()==1) {
			to_genotype = combined0_list.get(0);
			}
		else
			{
			to_genotype = invoke_combine(i,combined0_list,this.bedPath);
			i++;
			for(Path p:combined0_list) {
				FilesDelete(p);
				FilesDelete(indexFor(p));
				}
			}
		invoke_genotype(to_genotype,bedPath,this.outputFile);
		FilesDelete(to_genotype);
		FilesDelete(indexFor(to_genotype));
		
		return 0;
		}
	catch(Throwable err) {
		err.printStackTrace();
		return -1;
		}
	}
public static void main(String[] args) {
	new GatkHaplotypeCaller().instanceMainWithExit(args);
	}
}
