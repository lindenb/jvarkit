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
package com.github.lindenb.jvarkit.tools.liftover;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
/** 
 
## Example

 ```
 $ wget -q  -O - "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" | gunzip -c | java -jar dist/convertliftoverchain.jar -R1 src/test/resources/human_b37.dict  | grep chain -A 2 -m 4

chain 20851231461 1 249250621 + 10000 249240621 1 248956422 + 10000 248946422 2
167376	50041	80290
40302	253649	288020
--
chain 51897113 1 249250621 + 144274511 149034137 1 248956422 - 99285234 105680422 98
365	72	72
494	23	23
--
chain 24611930 1 249250621 + 206072707 206332221 1 248956422 - 42687778 42947276 254
33374	0	21
6406	0	1
--
chain 14571217 1 249250621 + 317719 471368 1 248956422 - 248454805 248608454 495
153649

[INFO][ConvertLiftOverChain]number of lines skipped 3346
[INFO][ConvertLiftOverChain]unmatched source contigs chr11_gl000202_random; chr17_ctg5_hap1; chr17_gl000203_random; chr17_gl000204_random; chr17_gl000205_random; chr17_gl000206_random; chr18_gl000207_random; chr19_gl000208_random; chr19_gl000209_random; chr1_gl000191_random; chr1_gl000192_random; chr21_gl000210_random; chr4_ctg9_hap1; chr4_gl000193_random; chr4_gl000194_random; chr6_apd_hap1; chr6_cox_hap2; chr6_dbb_hap3; chr6_mann_hap4; chr6_mcf_hap5; chr6_qbl_hap6; chr6_ssto_hap7; chr7_gl000195_random; chr8_gl000196_random; chr8_gl000197_random; chr9_gl000198_random; chr9_gl000199_random; chr9_gl000200_random; chr9_gl000201_random; chrUn_gl000211; chrUn_gl000212; chrUn_gl000213; chrUn_gl000214; chrUn_gl000215; chrUn_gl000216; chrUn_gl000217; chrUn_gl000218; chrUn_gl000219; chrUn_gl000220; chrUn_gl000221; chrUn_gl000222; chrUn_gl000223; chrUn_gl000224; chrUn_gl000225; chrUn_gl000226; chrUn_gl000227; chrUn_gl000228; chrUn_gl000229; chrUn_gl000230; chrUn_gl000231; chrUn_gl000232; chrUn_gl000233; chrUn_gl000234; chrUn_gl000235; chrUn_gl000236; chrUn_gl000237; chrUn_gl000238; chrUn_gl000239; chrUn_gl000240; chrUn_gl000241; chrUn_gl000242; chrUn_gl000243; chrUn_gl000244; chrUn_gl000245; chrUn_gl000246; chrUn_gl000247; chrUn_gl000248; chrUn_gl000249
[INFO][ConvertLiftOverChain]unmatched dest contigs chr11_KI270927v1_alt; chr12_GL877875v1_alt; chr14_GL000009v2_random; chr14_KI270726v1_random; chr14_KI270846v1_alt; chr15_KI270849v1_alt; chr15_KI270850v1_alt; chr15_KI270851v1_alt; chr15_KI270852v1_alt; chr17_KI270857v1_alt; chr17_KI270860v1_alt; chr17_KI270862v1_alt; chr17_KI270909v1_alt; chr19_KI270938v1_alt; chr1_KI270706v1_random; chr1_KI270711v1_random; chr1_KI270712v1_random; chr1_KI270765v1_alt; chr1_KI270766v1_alt; chr22_KI270731v1_random; chr22_KI270875v1_alt; chr22_KI270879v1_alt; chr22_KI270928v1_alt; chr2_GL383522v1_alt; chr2_KI270772v1_alt; chr2_KI270773v1_alt; chr2_KI270776v1_alt; chr2_KI270894v1_alt; chr4_GL000008v2_random; chr6_KI270801v1_alt; chr7_KI270803v1_alt; chr8_KI270811v1_alt; chr8_KI270821v1_alt; chr9_KI270719v1_random; chrUn_KI270742v1
```

 */
@Program(
		name="convertliftoverchain",
		description="Convert the contigs in a liftover chain to match another REFerence. (eg. to remove chr prefix, unknown chromosomes etc...)",
		keywords={"chain","liftover"},
		creationDate="20190409",
		modificationDate="20190409"
		)
public class ConvertLiftOverChain extends Launcher {
	private static final Logger LOG = Logger.build(ConvertLiftOverChain.class).make();

	@Parameter(names={"-o","--output"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R1","--ref1"},description="Source chain REFference mapping. Default : no conversion. "+ContigNameConverter.OPT_DICT_OR_MAPPING_FILE_DESC)
	private Path refFile1 = null;
	@Parameter(names={"-R2","--ref2"},description="Destination chain REFference mapping. Default : no conversion. "+INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path refFile2 = null;
	
	@Override
	public int doWork(final List<String> args) {
		PrintWriter out = null;
		BufferedReader in = null;
		final Pattern splitter = Pattern.compile("\\s");
		try {			
			final ContigNameConverter convert1 = this.refFile1==null ? ContigNameConverter.getIdentity() : ContigNameConverter.fromPathOrOneDictionary(this.refFile1);
			final ContigNameConverter convert2 = this.refFile2==null ? ContigNameConverter.getIdentity() : ContigNameConverter.fromPathOrOneDictionary(this.refFile2);
			
			final Set<String> notFound1 = new TreeSet<>();
			final Set<String> notFound2 = new TreeSet<>();
			
			in=super.openBufferedReader(oneFileOrNull(args));
			out = super.openPathOrStdoutAsPrintWriter(this.outputFile);
			int num_ignored=0;
			boolean valid=false;
			String line;
			while((line=in.readLine())!=null)
				{
				if(line.startsWith("chain")) {
					valid = false;
					final String chainFields[]=splitter.split(line);
					 if (chainFields.length != 13) {
						 out.close();out=null;
						 in.close();in=null;
						 throw new JvarkitException.TokenErrors(13, chainFields);
					 	}
					 final String fromSequenceName = chainFields[2];
					 String ctg = convert1.apply(fromSequenceName);
					 if(StringUtils.isBlank(ctg)) {
						 notFound1.add(fromSequenceName);
						 ++num_ignored;
						 continue;
					 	}
					 chainFields[2] = ctg;
					 final String toSequenceName = chainFields[7];
					 ctg = convert2.apply(toSequenceName);
					 if(StringUtils.isBlank(ctg)) {
						 notFound2.add(toSequenceName);
						 ++num_ignored;
						 continue;
					 	}
					 chainFields[7] = ctg;
					valid=true;
					out.println(String.join(" ", chainFields));
					}
				else if(valid)
					{
					out.println(line);
					}
				else
					{
					++num_ignored;
					}
				}
			out.flush();
			out.close();
			out= null;
			in.close();
			in = null;
			LOG.info("number of lines skipped "+num_ignored);
			LOG.info("unmatched source contigs "+notFound1.stream().collect(Collectors.joining("; ")));
			LOG.info("unmatched dest contigs "+notFound2.stream().collect(Collectors.joining("; ")));
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			CloserUtil.close(in);
			}
		}
	
	public static void main(final String[] args)
		{
		new ConvertLiftOverChain().instanceMainWithExit(args);
		}
}
