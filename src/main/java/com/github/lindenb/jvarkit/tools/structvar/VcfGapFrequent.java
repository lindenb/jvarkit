package com.github.lindenb.jvarkit.tools.structvar;

import java.io.Closeable;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.log.ProgressFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory;
import com.github.lindenb.jvarkit.util.vcf.AFExtractorFactory.AFExtractor;
import htsjdk.variant.vcf.VCFIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
BEGIN_DOC

## Example

```
$ find /commun/data/pubdb/broadinstitute.org/gnomad/release-170228/ -name "*.vcf.gz" > gnomad.list

$ java -jar dist/vcfgapfrequent.jar --fields "AF_POPMAX" -C -D gnomad.list genotyped.22.vcf.gz  > out.bed

```

END_DOC
 */

@Program(name="vcfgapfrequent",
description="Filter VCF annotated with external (AF or AC/AN) frequency information like vcfgnomad",
keywords={"vcf"}
)

public class VcfGapFrequent extends Launcher {
	private static final Logger LOG = Logger.build(VcfGapFrequent.class).make();
	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private File outputFile = null;

	@Parameter(names={"-D","--database"},description=
			"VCF database(s) used as a reference of frequent variants. "
			+ "One file ending with '.list' is interpretted as a list of path.")
	private List<String> databases = new ArrayList<>();

	@Parameter(names={"-f","--fields"},description="Where to peek the frequencies from the database."+AFExtractorFactory.OPT_DESC)
	private String user_fields_str = "AC/AN;AF";
	@Parameter(names={"-C","-skip","--skip"},description="Skip missing whole chromosome")
	private boolean skip_missing_whole_contig = false;
	@Parameter(names={"-m","--min-length"},description="Min segment length of NO_CALL")
	private int min_no_call_length = 1_000;
	@Parameter(names={"-M","--max-length"},description="Max segment length of NO_CALL")
	private int max_no_call_length = 1_000_000;
	@Parameter(names={"-t","--max","--treshold"},description=
			"Allele Frequency Treshold. Only Variant in database(s) with extracted AF greater than this value are considered.")
	private double af_treshold = 0.4;

	
	private PrintWriter out;
	private final List<FreqentDB> freqentDBs = new ArrayList<>();
	private SAMSequenceDictionary dict = null;
	
	private class FreqentDB
		implements Closeable
		{
		final VCFFileReader vcfFileReader;
		final ContigNameConverter converter;
		final List<AFExtractor> afExtractors = new ArrayList<>();

		
		
		FreqentDB(final File f) {
			this.vcfFileReader = new VCFFileReader(f,true);
			final VCFHeader header= this.vcfFileReader.getFileHeader();
			this.afExtractors.addAll(new AFExtractorFactory().parseFieldExtractors(user_fields_str));
			
			this.afExtractors.removeIf(EX->{
				if(!EX.validateHeader(header))
					{
					LOG.error("ignoreing "+EX+" for "+f+" because it's not valid");
					return true;
					}
				return false;
				});
			
			final SAMSequenceDictionary d = header.getSequenceDictionary();
			if(d==null)
				{
				converter = ContigNameConverter.getIdentity();
				}
			else
				{
				converter = ContigNameConverter.fromOneDictionary(d);
				}
			}
		@Override
		public void close() {
			CloserUtil.close(this.vcfFileReader);
			}
		}
	
	
	
	
	private class Handler
		{
		/** if null: no genotype is considered (whole VCF)*/
		final String sampleName;
		int prev_tid = 0;
		int prev_pos1 = 1;
		Handler(final String sampleName) {
			this.sampleName = sampleName;
			}
		void scan(int tid,int start_incl_1,int end_exclusive)
			{
			final int length = (end_exclusive-start_incl_1);
			if(length<=1) return;
			if(length<=min_no_call_length) return;
			if(length>=max_no_call_length) return;
						
			final SAMSequenceRecord ssr = dict.getSequence(tid);
			
			if(start_incl_1==1 && end_exclusive==ssr.getSequenceLength() && skip_missing_whole_contig)
				{
				return;
				}
			//LOG.debug("Scanning "+ssr.getSequenceName()+":"+start_incl_1+"-"+end_exclusive+" for "+this.sampleName);
			final Set<Integer> found_pos = new HashSet<>();
			for(final FreqentDB db:freqentDBs)
				{
				if(db.afExtractors.isEmpty()) continue;
				final String newCtg = db.converter.apply(ssr.getSequenceName());
				if(StringUtil.isBlank(newCtg)) continue;
				final CloseableIterator<VariantContext> iter = db.vcfFileReader.query(newCtg, start_incl_1,end_exclusive-1);
				while(iter.hasNext())
					{
					final VariantContext ctx = iter.next();
					if(found_pos.contains(ctx.getStart())) continue;
					final OptionalDouble maxAf= db.afExtractors.
							stream().flatMap(EX->EX.parse(ctx).stream()).
							mapToDouble(AF->AF==null ?0:AF.doubleValue()).
							filter(V-> V>= af_treshold).
							max();//pas noneMatch please
					if(!maxAf.isPresent()) continue;
					found_pos.add(ctx.getStart());
					}
				iter.close();
				}
			if(found_pos.isEmpty()) return;
			out.print(ssr.getSequenceName());
			out.print("\t");
			out.print(start_incl_1-1);
			out.print("\t");
			out.print(end_exclusive);
			out.print("\t");
			out.print(length);
			out.print("\t");
			out.print(this.sampleName==null?"*":this.sampleName);
			out.print("\t");
			out.print(found_pos.size());
			out.print("\t");
			out.print((double)found_pos.size()/(double)length);
			out.println();
			}
		
		void visit(final VariantContext ctx) {
			if(this.sampleName!=null)
				{
				final Genotype gt = ctx.getGenotype(this.sampleName);
				if(!(gt.isNoCall() || gt.getAlleles().stream().noneMatch(A->A.equals(Allele.SPAN_DEL)))) return;
				}
			final int ctx_tid = dict.getSequenceIndex(ctx.getContig());
			if(ctx_tid==-1) throw new JvarkitException.ContigNotFoundInDictionary(ctx.getContig(), dict);
			if(this.prev_tid> ctx_tid) {
				throw new IllegalStateException("Bad sort order "+ctx);
				}
			// finish current segment
			if(this.prev_tid<ctx_tid)
				{
				scan(this.prev_tid,this.prev_pos1,dict.getSequence(this.prev_tid).getSequenceLength());
				//go to next segment
				this.prev_pos1 = 1;
				this.prev_tid++;
				}
			
			// scan whole missing contig
			while(this.prev_tid<ctx_tid)
				{
				scan(this.prev_tid,1,dict.getSequence(this.prev_tid).getSequenceLength());
				this.prev_tid++;
				this.prev_pos1 = 1;
				}
			
			scan(ctx_tid,this.prev_pos1,ctx.getStart());
			this.prev_tid = ctx_tid;
			this.prev_pos1 = ctx.getEnd();
			}
		
		void finish() {
			scan(this.prev_tid,this.prev_pos1,dict.getSequence(this.prev_tid).getSequenceLength());
			this.prev_tid++;
			
			//scan remaining whole contigs
			while(this.prev_tid< dict.size())
				{
				scan(this.prev_tid,1,dict.getSequence(this.prev_tid).getSequenceLength());
				this.prev_tid++;
				this.prev_pos1=1;
				}
			}
		}
	

	
	@Override
	public int doWork(final List<String> args) {
		if(af_treshold<0 || this.af_treshold>1) 
			{
			LOG.error("bad treshold "+this.af_treshold);
			return -1;
			}
		
		VCFIterator in = null;
		try 
			{
			
			this.freqentDBs.addAll( IOUtils.unrollFiles2018(this.databases).
					stream().
					map(F->new FreqentDB(F)).collect(Collectors.toList())
					);
			if(this.freqentDBs.isEmpty())
				{
				LOG.error("no database was defined");
				return -1;
				}
			
			if(this.freqentDBs.stream().flatMap(DB->DB.afExtractors.stream()).count()==0L)
				{
				LOG.error("No AF extractor is suitable for any databases.");
				return -1;
				}
			
			in= super.openVCFIterator(oneFileOrNull(args));
			final VCFHeader header = in.getHeader();
			
			
			this.dict = header.getSequenceDictionary();
			if(this.dict==null || this.dict.isEmpty())
				{
				LOG.error(JvarkitException.VcfDictionaryMissing.getMessage("input"));
				return -1;
				}
			this.out= super.openFileOrStdoutAsPrintWriter(this.outputFile);
			
			final List<Handler> handlers = new ArrayList<>(header.getNGenotypeSamples()+1);
			handlers.add(new Handler(null));
			header.getSampleNamesInOrder().stream().map(S->new Handler(S)).forEach(H->handlers.add(H));
			
			final ProgressFactory.Watcher<VariantContext> progress = ProgressFactory.newInstance().
					dictionary(header).
					validatingSortOrder(true).
					logger(LOG).
					build();
			
			while(in.hasNext())
				{
				final VariantContext ctx =progress.apply(in.next());
				handlers.stream().forEach(H->H.visit(ctx));
				}
			handlers.stream().forEach(H->H.finish());
			progress.close();
			
			
			out.flush();
			out.close();
			out=null;
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(in);
			CloserUtil.close(out);
			CloserUtil.close(this.freqentDBs);
			}
		}
	
	public static void main(final String[] args) {
		new VcfGapFrequent().instanceMainWithExit(args);
	}
	
}
