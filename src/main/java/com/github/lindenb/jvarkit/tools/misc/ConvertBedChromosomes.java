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


History:
* 2014 creation

*/
package com.github.lindenb.jvarkit.tools.misc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import htsjdk.samtools.util.CloserUtil;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.NoSplitter;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;

/**
BEGIN_DOC


## Example

```bash
$   curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz |\
    gunzip -c | \
    java -jar dist/bedrenamechr.jar -f src/main/resources/chromnames/hg19_to_g1kv37.tsv -c 2 |\
   tail


uc011nca.2	Y	+	59213948	59276439	59230880	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222281,59230919,59233257,59252550,59272463,59276439,	P51809	uc011nca.2
uc004fxl.3	Y	+	59213948	59276439	59222135	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222216,59230919,59233257,59252550,59272463,59276439,	P51809-3	uc004fxl.3
uc004fxk.3	Y	+	59213948	59276439	59222135	59274809	7	59213948,59222126,59228291,59230781,59233166,59272370,59274552,	59214117,59222281,59228349,59230919,59233257,59272463,59276439,	P51809-2	uc004fxk.3
uc011ncb.2	Y	+	59213948	59276439	59222274	59274621	7	59213948,59222126,59230781,59233166,59252482,59272370,59274552,	59214117,59222277,59230919,59233257,59252550,59272463,59276439,	B4DE96	uc011ncb.2
uc010nxr.2	Y	+	59330251	59340490	59335611	59340461	9	59330251,59334000,59335576,59336119,59336347,59337090,59337948,59338753,59340193,	59330458,59334179,59335690,59336231,59336526,59337236,59338150,59338859,59340490,	B4E011	uc010nxr.2
uc004fxm.1	Y	+	59330251	59343488	59330414	59342523	9	59330251,59334078,59335552,59336119,59336354,59337119,59337948,59338753,59342486,	59330458,59334179,59335690,59336231,59336526,59337236,59338150,59338859,59343488,	B9ZVT0	uc004fxm.1
uc004fxn.1	Y	+	59330251	59343488	59330430	59343080	9	59330251,59335576,59336119,59336347,59337090,59337948,59338753,59340193,59342486,	59330458,59335690,59336231,59336526,59337236,59338150,59338859,59340278,59343488,	Q01113	uc004fxn.1
uc004fxo.1	Y	+	59352972	59356131	59353631	59356130	7	59352972,59354351,59354669,59354993,59355369,59355682,59355972,	59353819,59354463,59354816,59355130,59355505,59355884,59356131,	I3L0A4	uc004fxo.1
uc022cpg.1	Y	+	59354984	59358336	59355427	59358045	7	59354984,59355369,59355682,59355972,59356790,59357702,59357911,	59355130,59355505,59355884,59356131,59356943,59357771,59358336,	Q9NQA3	uc022cpg.1
uc011ncc.1	Y	-	59358328	59360854	59358328	59358328	3	59358328,59360006,59360500,	59359508,59360115,59360854,	uc011ncc.1

```

END_DOC

 */
@Program(
		name="bedrenamechr",
		description="Convert the names of the chromosomes in a Bed file",
		keywords={"bed","chromosome","contig","convert"},
		modificationDate="20190503"
		)
public class ConvertBedChromosomes
	extends Launcher
	{
	private static final Logger LOG = Logger.build(ConvertBedChromosomes.class).make();
	
	private  enum OnNotFound{RAISE_EXCEPTION,SKIP,RETURN_ORIGINAL};

	
	@Parameter(names={"-convert","--convert"},description="What should I do when  a converstion is not found")
	private OnNotFound onNotFound=OnNotFound.RAISE_EXCEPTION;
	@Parameter(names={"-f","--mapping","-m","-R"},description="load a custom name mapping."+ ContigNameConverter.OPT_DICT_OR_MAPPING_FILE_DESC,required=true)
	private Path mappingFile=null;
	@Parameter(names={"-c","--column"},description="1-based chromosome column(s), multiple separated by commas",splitter=NoSplitter.class)
	private String chromColumnsStr="1";
	@Parameter(names={"-o","--out"},description=OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile= null;
	@Parameter(names={"-s","--header"},description="Ignore lines starting with this java regular expression")
	private String ignoreLinesPattern="(#|browser|track)";
	@Parameter(names={"-d","--delim"},description="field delimiter.")
	private String delimStr="\t";

	
	private ContigNameConverter customMapping=ContigNameConverter.getIdentity();
	private Set<String> unmappedChromosomes=new HashSet<String>();
	private int chromColumns0[];
	
	
	@SuppressWarnings("resource")
	protected int doWork(final BufferedReader in,PrintStream out)
			throws IOException
		{
		String line;
		final  CharSplitter tab = CharSplitter.of(this.delimStr);
		final Pattern headerRegex = Pattern.compile(this.ignoreLinesPattern);
		
		while((line=in.readLine())!=null)
			{
			final Matcher match =  headerRegex.matcher(line);
			if(match.find() && match.start()==0)
				{
				out.println(line);
				continue;
				}
			boolean ok=true;
			final String tokens[]=tab.split(line);
			for(final int chromColumn0: this.chromColumns0) 
				{
				if(chromColumn0<0 || chromColumn0 >=tokens.length) {
					ok=false;
					continue;
					}
				final String chrom= this.customMapping.apply(tokens[chromColumn0]);
				if(StringUtils.isBlank(chrom)) {
					if(this.unmappedChromosomes.add(tokens[chromColumn0])) {
						LOG.warn("cannot convert "+ tokens[chromColumn0]);
						}
					ok=false;
					continue;
					}
				tokens[chromColumn0]=chrom;
				}
			if(!ok) {
				switch(this.onNotFound) {
					case RAISE_EXCEPTION:  throw new IOException("Cannot convert chrom in : "+line);
					case SKIP: continue;
					case RETURN_ORIGINAL: break;
					}
				}
			out.println(String.join(String.valueOf(tab.getDelimiter()), tokens));
			}
		out.flush();
		return 0;
		}
	
	@Override
	public int doWork(final List<String> args) {
		
		PrintStream out=null;
		try
			{
			this.customMapping=ContigNameConverter.fromPathOrOneDictionary(this.mappingFile);
			
			this.chromColumns0 = Arrays.stream(CharSplitter.COMMA.split(this.chromColumnsStr)).
				filter(S->!StringUtils.isBlank(S)).
				map(S->Integer.parseInt(S)-1)./* convert to 0-based */
				collect(Collectors.toCollection(TreeSet::new)).
				stream().
				mapToInt(I->I.intValue()).
				toArray();
			
			if(this.chromColumns0.length==0) {
				LOG.error("No column defined");
				return -1;
				}
			
			out = super.openPathOrStdoutAsPrintStream(this.outputFile);
			if(args.isEmpty())
				{
				doWork(IOUtils.openStreamForBufferedReader(stdin()), out);
				}
			else
				{
				for(final String filename:args) {
					try(final BufferedReader br=IOUtils.openURIForBufferedReading(filename)) {
						doWork(br, out);
						}
					}
				}
			if(!unmappedChromosomes.isEmpty())
				{
				LOG.warning("Unmapped chromosomes:"+unmappedChromosomes);
				}
			out.flush();
			return 0;
			}
		catch(final Throwable err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			CloserUtil.close(out);
			}
		}
	

	public static void main(final String[] args)
		{
		new ConvertBedChromosomes().instanceMainWithExit(args);
		}
	}
