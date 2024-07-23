package com.github.lindenb.jvarkit.tools.gatk;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.gatk.Gatk4Proxy;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.vcf.VCFUtils;

import htsjdk.samtools.util.FileExtensions;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFIterator;
import htsjdk.variant.vcf.VCFIteratorBuilder;

public class GatkHaplotypeCaller extends Launcher {
	private static final Logger LOG = Logger.build(GatkHaplotypeCaller.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description=INDEXED_FASTA_REFERENCE_DESCRIPTION,required=true)
	private Path reference = null;
	@Parameter(names={"-bed","--bed"},description="restrict to bed",required=true)
	private Path bedPath = null;
	@Parameter(names={"-dbsnp","--dbsnp"},description="path to dbsnp")
	private Path dbsnp = null;

	private int mapq=10;
	
	private final Gatk4Proxy gatkEngine = Gatk4Proxy.getInstance().orElse(null);

	private void execute(final List<String> argv) throws Exception {
		this.gatkEngine.execute(argv);
		}

	private void fill_cmd(List<String> cmd,final Path tmpDir,final Path bed) {
		cmd.add("-R");
		cmd.add(this.reference.toString());
		//cmd.add("--verbosity");
		//cmd.add("ERROR");
		if(bed!=null) {
			cmd.add("-L");
			cmd.add(bed.toString());
			}
		cmd.add("--tmp-dir");
		cmd.add(tmpDir.toString());
		if(this.dbsnp!=null) {
			cmd.add("--dbsnp");
			cmd.add(String.valueOf(this.dbsnp));
			}
		}

	private Path invoke_gvcf(final Path tmpDir,int idx,final Path bam,final Path bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("hc.%03d.g" + FileExtensions.COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("HaplotypeCaller");
		fill_cmd(cmd,tmpDir,bed);
		cmd.add("-I");
		cmd.add(bam.toString());
		cmd.add("-O");
		cmd.add(outPath.toString());
		cmd.add("-ERC");
		cmd.add("GVCF");
		if(this.mapq>0) {
			cmd.add("--minimum-mapping-quality");
			cmd.add(String.valueOf(this.mapq));
			}
		
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-G");cmd.add("StandardHCAnnotation");
		
		execute(cmd);
		return outPath;
		}

	private Path invoke_combine(final Path tmpDir,int idx,final List<Path> gvcfs,final Path bed) throws Exception {
		final Path outPath = tmpDir.resolve(String.format("combine.%03d.g" + FileExtensions.COMPRESSED_VCF, idx));
		final List<String> cmd= new ArrayList<>();
		cmd.add("CombineGVCFs");
		fill_cmd(cmd,tmpDir,bed);
		for(Path gvcf:gvcfs) {
			cmd.add("-V");
			cmd.add(gvcf.toString());
			}
		cmd.add("-O");
		cmd.add(outPath.toString());
		
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");

		
		execute(cmd);
		return outPath;
		}
	
	private void invoke_genotype(final Path tmpDir,final Path gvcf,final Path bed,final Path outPath) throws Exception {
		final List<String> cmd= new ArrayList<>();
		final Path vcfOut;
		if(outPath==null) {
			vcfOut = Files.createTempFile(tmpDir,"hc.",".vcf.gz");
			}
		else
			{
			vcfOut=outPath;
			}
		cmd.add("GenotypeGVCFs");
		fill_cmd(cmd,tmpDir,bed);
		cmd.add("-V");
		cmd.add(gvcf.toString());
		cmd.add("-G");cmd.add("StandardAnnotation");
		cmd.add("-G");cmd.add("AS_StandardAnnotation");
		cmd.add("-O");
		cmd.add(vcfOut.toString());
		execute(cmd);
		if(outPath==null) {
			try(VCFIterator r=new VCFIteratorBuilder().open(vcfOut)) {
				try(VariantContextWriter w= VCFUtils.createVariantContextWriterToOutputStream(stdout())) {
					w.writeHeader(r.getHeader());
					VCFUtils.copyVariantsTo(r, w);
					}
				}
			}
		}


	
@Override
public int doWork(List<String> args) {
	try {
		if(gatkEngine==null) {
			System.err.println("gatk4 engine is not compiled");
			return -1;
			}
		final List<Path> bams = IOUtils.unrollPaths(args);
		
		if (bams.isEmpty()) {
        	LOG.error("No Bam was provided");
        	return -1;
        	}
		
		
		Path tmpDir = Files.createTempDirectory("tmp");
		final List<Path> gvcfs_list = new ArrayList<>(bams.size());
		for(int i=0;i< bams.size();i++) {
			LOG.info("("+(i+1)+"/"+bams.size()+" "+bams.get(i));
			final Path g_vcf_gz = invoke_gvcf(tmpDir,i,bams.get(i),this.bedPath);
			gvcfs_list.add(g_vcf_gz);
			}
			
		final int sqrt = Math.max(20, (int)Math.sqrt(gvcfs_list.size()));
		
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
				final Path combined = invoke_combine(tmpDir,i,L, this.bedPath);
				combined0_list.add(combined);
				}
			i++;
			}
		Path to_genotype;
		if(combined0_list.size()==1) {
			to_genotype = combined0_list.get(0);
			}
		else
			{
			to_genotype = invoke_combine(tmpDir,i,combined0_list,this.bedPath);
			i++;
			}
		invoke_genotype(tmpDir,to_genotype,bedPath,this.outputFile);

		
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
