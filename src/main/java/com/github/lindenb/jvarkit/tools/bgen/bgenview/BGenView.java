package com.github.lindenb.jvarkit.tools.bgen.bgenview;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.function.UnaryOperator;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bgen.BGenGenotype;
import com.github.lindenb.jvarkit.bgen.BGenReader;
import com.github.lindenb.jvarkit.bgen.BGenUtils;
import com.github.lindenb.jvarkit.bgen.BGenVariant;
import com.github.lindenb.jvarkit.bgen.BGenWriter;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.interval.TargetsParameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Md5CalculatingOutputStream;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFStandardHeaderLines;


/**
BEGIN_DOC

## Example

```

```

END_DOC
 */
@Program(name="bgenview",
description="utilities for bgen",
keywords={"bgen","vcf"},
creationDate="20250509",
modificationDate="20250509",
jvarkit_amalgamion =  true,
generate_doc = false
)
public class BGenView extends Launcher {
	private static final Logger LOG = Logger.of(BGenView.class).setDebug();


	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	
	@Parameter(names={"--compression"},description="Compression type")
	private BGenUtils.Compression compression = BGenUtils.Compression.ZLIB;
	
	//@Parameter(names={"--layout"},description="BGen Layout")
	//private BGenUtils.Layout layout = BGenUtils.Layout.e_Layout2;
	
	@Parameter(names={"--dict"},description="Use a dictionary to rename the contig. If a contig is missing, just skip the variant. "+ DICTIONARY_SOURCE)
	private Path pathDict=null;

	@Parameter(names={"--anonymous"},description="make anonymous, hide sample names")
	private boolean anonymous_samples=false;
	@Parameter(names={"--bits"},description="number of bits per prob. <=0 use incoming variant capacity")
	private int nbits=-1;

	@ParametersDelegate
	TargetsParameter regionFilesParameter = new TargetsParameter();

	
	public void writeVariant(
			final BGenVariant source_ctx,
			String contig,int pos,String variantId,String rsId,List<String> alleles,
			final int[] sample_new_index,
			BGenWriter w
			) throws IOException{
		
		if(variantId.isEmpty()) {
			final Md5CalculatingOutputStream md5 = new Md5CalculatingOutputStream(OutputStream.nullOutputStream(), (Path)null);
			md5.write(contig.getBytes(BGenUtils.ENCODING));
			md5.write(String.valueOf(pos).getBytes(BGenUtils.ENCODING));
			for(String a:alleles) md5.write(a.getBytes(BGenUtils.ENCODING));
			md5.close();
			variantId = md5.md5().substring(0,8);
			}
		
		w.writeVariant(contig,pos,variantId,rsId,alleles);
		w.setPhased(source_ctx.isPhased());
		for(int x=0;x < source_ctx.getNGenotypes();++x)  {
			if(sample_new_index[x]==-1) continue;//this sample is removed
			final BGenGenotype gt = source_ctx.getGenotype(x);
			w.setGenotype(
					sample_new_index[x],
					gt.getPloidy(),
					gt.isMissing(),
					gt.getProbs()
					);
			}
		w.writeGenotypes();
		}
	
	private void writeVariant(
			UnaryOperator<String> ctgRename,
			BGenVariant ctx,
			final int[] sample_new_index,
			BGenWriter w
			) throws IOException {
		final String contig = ctgRename.apply(ctx.getContig());
		if(StringUtils.isBlank(contig)) return;
		writeVariant(
			ctx,
			contig, ctx.getPosition(), ctx.getId(),ctx.getRsId(),ctx.getAlleles(),
			sample_new_index,
			w);
		
		}
	
	@Override
	public int doWork(final List<String> args) {
		try {
			final SAMSequenceDictionary dict=(this.pathDict==null?null:new SequenceDictionaryExtractor().extractRequiredDictionary(this.pathDict));
			final UnaryOperator<String> renameCtg =  (dict==null?null:ContigNameConverter.fromOneDictionary(dict));
			final Predicate<Locatable> accept = regionFilesParameter.makePredicate();
			final String input = oneFileOrNull(args);
			try(BGenReader r=(input==null || input.equals("-")?
					new BGenReader(stdin()):
					new BGenReader(input))) {
					try(BGenWriter w = (this.outputFile==null?
							new BGenWriter(stdout()):
							new BGenWriter(this.outputFile)
							)) {
					w.setCompression(this.compression);
					w.setAnonymousSamples(this.anonymous_samples);
					w.setLayout(r.getHeader().getLayout());
					w.setNBits(this.nbits);
					
					final List<String> input_samples= r.getHeader().getSamples();
					final List<String> output_samples= new ArrayList<>(input_samples.size());

					final int[] sample_new_index = new int[input_samples.size()];
					
					int j=0;
					for(int i=0;i< input_samples.size();i++) {
						//TODO add here sample filtering
						sample_new_index[i]=j;
						output_samples.add(input_samples.get(i));
						j++;
						}
					
					
					
					if(r.getHeader().hasAnonymousSamples()) {
						w.writeHeader(output_samples.size());
						}
					else
						{
						w.writeHeader(output_samples);
						}
					
					for(;;) {
						BGenVariant ctx = r.readVariant();
						if(ctx==null) break;
						if(!accept.test(ctx)) {
							r.skipGenotypes();
							continue;
							}
						ctx = r.readGenotypes();
						
						if(LOG.isDebug()) {
							LOG.debug(ctx);
						}
						
						w.setNBits(this.nbits<=0?ctx.getBitsPerProb():this.nbits);
						writeVariant(renameCtg,ctx, sample_new_index, w);
						}
					
					}
				}
			return 0;
			}
		catch(Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally 
			{
			
			}
		}
	
	public static void main(String[] args) {
		new BGenView().instanceMainWithExit(args);
	}

}
