package com.github.lindenb.jvarkit.tools.bgen.bgen2vcf;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.bgen.BGenGenotype;
import com.github.lindenb.jvarkit.bgen.BGenReader;
import com.github.lindenb.jvarkit.bgen.BGenUtils;
import com.github.lindenb.jvarkit.bgen.BGenVariant;
import com.github.lindenb.jvarkit.dict.SequenceDictionaryExtractor;
import com.github.lindenb.jvarkit.interval.TargetsParameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.tools.bgen.bgenview.BGenView;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.variant.variantcontext.writer.WritingVariantsDelegate;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
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
@Program(name="bgen2vcf",
description="convert bgen to vcf",
keywords={"bgen","vcf"},
creationDate="20250508",
modificationDate="20250509",
jvarkit_amalgamion =  true,
generate_doc = false
)
public class BGenToVcf extends Launcher {

	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile=null;
	@Parameter(names={"-R","--reference"},description=DICTIONARY_SOURCE)
	private Path dictPath=null;
	@Parameter(names={"-G","--no-genotype"},description="skip genotypes")
	private boolean skip_genotypes=false;
	@ParametersDelegate
	TargetsParameter regionFilesParameter = new TargetsParameter();
	@Parameter(names={"-N","--head"},description="limit to 'N' variant; negative==show all ")
	private long limit_n_variants=-1L;

	@ParametersDelegate
	private WritingVariantsDelegate writingVariantsDelegate=new WritingVariantsDelegate();

	
	private static final Logger LOG = Logger.build(BGenView.class).make();
	
	
	@Override
	public int doWork(List<String> args) {
		try {
			SAMSequenceDictionary dict = null;
			if(this.dictPath!=null) {
				dict = new SequenceDictionaryExtractor().extractRequiredDictionary(this.dictPath);
				}
			final Predicate<Locatable> accept = regionFilesParameter.
					setDictionary(dict).
					makePredicate();
			
			
			final String input = oneFileOrNull(args);
			try(BGenReader r=(input==null || input.equals("-")?
					new BGenReader(stdin()):
					new BGenReader(input))) {

				final Set<VCFHeaderLine> metaData = new HashSet<>();
				VCFStandardHeaderLines.addStandardFormatLines(metaData, true, VCFConstants.GENOTYPE_KEY);
				VCFStandardHeaderLines.addStandardInfoLines(metaData, true,VCFConstants.END_KEY);

				final VCFInfoHeaderLine BITS_info_header_line = new VCFInfoHeaderLine(
						"BITS",1,VCFHeaderLineType.Integer,"Number of bits used for storage in the bgen file");
				
				final VCFInfoHeaderLine ID_info_header_line = new VCFInfoHeaderLine(
						"ID2",1,VCFHeaderLineType.String,"bgen variant id");

				
				metaData.add(BGenUtils.GP_format_header_line);
				metaData.add(BGenUtils.HP_format_header_line);
				metaData.add(BITS_info_header_line);
				metaData.add(BGenUtils.OFFSET_format_header_line);
				
				metaData.add(new VCFHeaderLine("bgen.compression", r.getHeader().getCompression().name()));
				metaData.add(new VCFHeaderLine("bgen.n-variants", String.valueOf( r.getHeader().getNVariants())));
				metaData.add(new VCFHeaderLine("bgen.n-samples", String.valueOf( r.getHeader().getNSamples())));
				metaData.add(new VCFHeaderLine("bgen.layout", r.getHeader().getLayout().name()));
				metaData.add(new VCFHeaderLine("bgen.snps-offset", String.valueOf( r.getSnpsOffset())));
				metaData.add(new VCFHeaderLine("bgen.anonymous", String.valueOf( r.getHeader().hasAnonymousSamples())));

				
				final VCFHeader header;
				
				if(this.skip_genotypes) {
					header=new VCFHeader(metaData);
					}
				else
					{
					header=new VCFHeader(metaData,r.getHeader().getSamples());
					}
				if(dict!=null) header.setSequenceDictionary(dict);
				JVarkitVersion.getInstance().addMetaData(this, header);
				long n_variants=0L;
				try(VariantContextWriter w=this.writingVariantsDelegate.dictionary(dict).open(this.outputFile)) {
					w.writeHeader(header);
					
					for(;;) {
						if(limit_n_variants>=0 && n_variants>=limit_n_variants) break;
						
						BGenVariant ctx = r.readVariant();
						if(ctx==null) break;
						if(!accept.test(ctx)) {
							r.skipGenotypes();
							continue;
							}
						if(dict!=null && dict.getSequence(ctx.getContig())==null) {
							throw new JvarkitException.ContigNotFoundInDictionary(ctx.getContig(), dict);
							}
						n_variants++;
						
						final List<Allele> alleles = new ArrayList<>(ctx.getNAlleles());
						for(int a=0;a< ctx.getNAlleles();++a) {
							alleles.add(Allele.create(ctx.getAllele(a), a==0));
							}
						final int endPos = ctx.getStart() + alleles.get(0).length()-1;
						final VariantContextBuilder vcb=new VariantContextBuilder(input, 
								ctx.getContig(), 
								ctx.getStart(), 
								endPos, 
								alleles
								);
						
						final long offset = ctx.getOffset();
						if(offset>0L) {
							vcb.attribute(BGenUtils.OFFSET_format_header_line.getID(),
								String.valueOf(offset)/* as a String because it can be a int63 number, not an int32 */
								);
							}
						
						final String id = ctx.getId();
						if(!StringUtils.isBlank(id) && !id.equals(".")) {
							vcb.attribute(ID_info_header_line.getID(),id);
							}
						if(!StringUtils.isBlank(ctx.getRsId())) {
							vcb.id(ctx.getRsId());
							}
						if(ctx.getStart() < endPos) {
							vcb.attribute(VCFConstants.END_KEY, endPos);
							}
						if(this.skip_genotypes) {
							r.skipGenotypes();
							}
						else
							{
							ctx = r.readGenotypes();
							vcb.attribute(BITS_info_header_line.getID(), ctx.getBitsPerProb());
							
							final List<Genotype> genotypes = new ArrayList<>(ctx.getNGenotypes());
							for(int i=0;i< ctx.getNGenotypes();i++) {
								final BGenGenotype gt  =ctx.getGenotype(i);
								if(gt.isMissing()) {
									genotypes.add(GenotypeBuilder.createMissing(gt.getSample(), gt.getPloidy()));
									}
								else
									{
									// build from missing otherwise "java.lang.IllegalStateException: GTs cannot be missing for some samples if they are available for others in the record"
									final GenotypeBuilder gb=new GenotypeBuilder(GenotypeBuilder.createMissing(gt.getSample(), gt.getPloidy()));
									gb.phased(gt.isPhased());
									gb.attribute(BGenUtils.GP_format_header_line.getID(), gt.getProbs());
									genotypes.add(gb.make());
									}
								}
							vcb.genotypes(genotypes);
							}
						
						
						w.add(vcb.make());
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
		new BGenToVcf().instanceMainWithExit(args);
	}

}
