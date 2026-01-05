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
package com.github.lindenb.jvarkit.tools.gtf;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.jcommander.Launcher;
import com.github.lindenb.jvarkit.jcommander.Program;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

/**

BEGIN_DOC

```
$ java -jar ${JVARKIT_DIST}/gtf2bed.jar --columns gtf.feature,gene_name,gene_biotype,gene_id ~/src/jvarkit/src/test/resources/Homo_sapiens.GRCh37.87.gtf.gz | head | column -t
#chrom  start      end        gtf.feature      gene_name  gene_biotype    gene_id
1       120454175  120459317  exon             NOTCH2     protein_coding  ENSG00000134250
1       120454175  120612240  gene             NOTCH2     protein_coding  ENSG00000134250
1       120454175  120457928  three_prime_utr  NOTCH2     protein_coding  ENSG00000134250
1       120454175  120612240  transcript       NOTCH2     protein_coding  ENSG00000134250
1       120457928  120457931  stop_codon       NOTCH2     protein_coding  ENSG00000134250
1       120457931  120459317  CDS              NOTCH2     protein_coding  ENSG00000134250
1       120460287  120460385  CDS              NOTCH2     protein_coding  ENSG00000134250
1       120460287  120460385  exon             NOTCH2     protein_coding  ENSG00000134250
1       120461028  120461176  CDS              NOTCH2     protein_coding  ENSG00000134250
```

END_DOC
*/
@Program(
		name="gtf2bed",
		description="Convert GTF/GFF3 to BED.",
		creationDate="20220629",
		modificationDate="20220630",
		keywords= {"gtf","gff","gff3","bed"},
		jvarkit_amalgamion =  true,
		menu="GTF/GFF Manipulation"
		)
public class GtfToBed
	extends Launcher {
	private static final Logger LOG = Logger.of(GtfToBed.class);
	
	private final static String MY_COLS="gtf.source,gtf.feature,gtf.score,gtf.strand,gtf.frame";

	private final static String GTF_COLS = "ccds_id,exon_id,exon_number,exon_version,gene_biotype,gene_id,gene_name,gene_source,"
		+ "gene_version,havana_transcript,havana_transcript_version,protein_id,protein_version,tag,"
		+ "transcript_biotype,transcript_id,transcript_name,transcript_source,transcript_version";
	private final static String GFF_COLS= 
		"Alias,biotype,ccdsid,constitutive,description,ensembl_end_phase,"
			+ "ensembl_phase,exon_id,external_name,gene_id,havana_transcript,"
			+ "havana_version,ID,logic_name,Name,Parent,protein_id,rank,tag,"
			+ "transcript_id,version"
			;
	private final String ALL_COLS = MY_COLS+","+GTF_COLS+","+GFF_COLS;

	
	@Parameter(names={"-o","--output"},description= OPT_OUPUT_FILE_OR_STDOUT)
	private Path outputFile = null;
	@Parameter(names={"-R","--reference"},description= INDEXED_FASTA_REFERENCE_DESCRIPTION)
	private Path faidx = null;

	@Parameter(names={"-c","--columns"},description= "comma separated columns to be displayed")
	private String columnStr=ALL_COLS;

	@Parameter(names={"--grep"},description= "Check some identifiers are found in a column. syntax: <COLUMN>:<FILE_CONTAINING_THE_IDENTIFIERS>")
	private List<String> grepStr = new ArrayList<>();


	private static class Grep {
		final String column;
		final Path path;
		final Set<String> all;
		final Set<String> remains;
		Grep(final String column,final Path path) {
			this.column = column;
			this.path = path;
			try {
				this.all= Files.lines(this.path).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet());
				this.remains = new HashSet<>(this.all);
				}
			catch(IOException err) {
				throw new RuntimeIOException(err);
				}
			}
		boolean visit(Map<String,String> map) {
			if(!map.containsKey(this.column)) return false;
			final String value = map.get(this.column);
			if(!this.all.contains(value)) return false;
			this.remains.remove(value);
			return true;
			}
		boolean isSuccess() {
			if(this.remains.isEmpty()) return true;
			System.err.print("in file "+path+" for column "+this.column+" the following items where not found:");
			for(String s:this.remains) System.err.print(" "+s);
			System.err.println();
			return false;
			}
		}

	
	private	final Set<String> columns = new LinkedHashSet<>();
	private enum Type{ undefined,gtf,gff3}
	private Type gtype = Type.undefined;

	private void put(final Map<String,String> map,final String key,final String value) {
		if(key.isEmpty()) throw new IllegalArgumentException("key is empty ");
		if(!this.columns.contains(key)) return;
		if(map.put(key, value)!=null) {
			throw new IllegalArgumentException("duplicate key "+key+" in "+map);
			}
		}

	private Map<String,String> attributes(final String s) {
		final Map<String,String> map= new HashMap<>();
		if(s.equals(".")) return map;
		int i=0;
		
		while(i< s.length()) {
			/* skip ws */
			while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
			final StringBuilder key=new StringBuilder();
			while(i< s.length()) {
				if(this.gtype.equals(Type.undefined)) {
					if(Character.isWhitespace(s.charAt(i))) {
						this.gtype = Type.gtf;
						i++;
						break;
						}
					else if(s.charAt(i)=='=') {
						this.gtype = Type.gff3;
						i++;
						break;
						}
					}
				if(this.gtype.equals(Type.gtf) && Character.isWhitespace(s.charAt(i))) {
					i++;
					break;
					}
				else if(this.gtype.equals(Type.gff3) && s.charAt(i)=='=') {
					i++;
					break;
					}
				key.append(s.charAt(i));
				i++;
				}
			if(this.gtype.equals(Type.undefined)) throw new IllegalArgumentException("undefined gtf type with " + s);
			/* skip ws */
			while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
			if(i>=s.length()) throw new IllegalArgumentException("END of string after "+s.substring(0,i));
			

			boolean got_quote=false;
			if(this.gtype.equals(Type.gtf) && s.charAt(i)=='\"') {
				got_quote=true;
				i++;
				}
			
			
			final StringBuilder value=new StringBuilder();
			while(i< s.length() ) {
				if(got_quote && this.gtype.equals(Type.gtf) && s.charAt(i)=='\"') break;
				else if(!got_quote && this.gtype.equals(Type.gtf) && s.charAt(i)==';') break;
				else if(!got_quote && this.gtype.equals(Type.gtf) && Character.isWhitespace(s.charAt(i))) break;
				else if(this.gtype.equals(Type.gff3) && s.charAt(i)==';') break;
				
				if(s.charAt(i)=='\\') {
					if(i+1>=s.length()) throw new IllegalArgumentException("unknown escape sequence after "+s.substring(0,i));
					i++;
					switch(s.charAt(i)) {
						case '"': value.append("\"");break;
						case '\'': value.append("\'");break;
						case '\\': value.append("\\");break;
						case 't': value.append("\t");break;
						case 'n': value.append("\n");break;
						case '/': value.append("/");break;// got the string 'kinase\/fructos" in a gtf...
						default: throw new IllegalArgumentException("unknown escape sequence after "+s.substring(0,i));
						}
					}
				else
					{
					value.append(s.charAt(i));
					}
				i++;
				}
			if(got_quote && this.gtype.equals(Type.gtf)) {
				if(i>=s.length()) throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
				if(this.gtype.equals(Type.gtf) &&  s.charAt(i)!='\"')  throw new IllegalArgumentException("expected quote after "+s.substring(0,i));
				i++;
				}
			this.put(map,key.toString(), value.toString());
			
			/* skip ws */
			while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
			if(i<s.length()) {
				if(s.charAt(i)!=';')  throw new IllegalArgumentException("expected semicolon after "+s.substring(0,i));
				i++;
				}
			/* skip ws */
			while(i< s.length() && Character.isWhitespace(s.charAt(i))) i++;
			}
		return map;
		}

	
	@Override
	public int doWork(final List<String> args)  {
		try {
			final List<Grep> greps = new ArrayList<>();

			final String[] fixed_cols_prefixes=new String[]{"gtf.","gff.","gff3."};
			final ContigNameConverter contigNameConverter =  
					this.faidx==null?
					ContigNameConverter.getIdentity():
					ContigNameConverter.fromOneDictionary(SequenceDictionaryUtils.extractRequired(this.faidx))
					;
			this.columns.clear();
			this.columns.addAll(Arrays.asList(this.columnStr.split("[, \t;\\|]+")));
			this.columns.remove("");

			for(final String str:this.grepStr) {
				final int colon = str.indexOf(":");
				if(colon==-1) {
					LOG.error("expected ':' in grep expression "+str);
					return -1;
					}
				final String col = str.substring(0,colon);
				if(!this.columns.contains(col)) {
					LOG.error("column "+col+" not used by user for grep expression "+str +" ("+String.join(",",columns)+")");
					return -1;
					}
				final Path p = Paths.get(str.substring(colon+1));
				IOUtil.assertFileIsReadable(p);
				greps.add(new Grep(col,p));
				}
			
			final String input = oneFileOrNull(args);
			try(InputStream in1 = input==null?stdin():Files.newInputStream(Paths.get(input))) {
				final InputStream in = IOUtils.uncompress(in1);
				final BufferedReader br = new BufferedReader(new InputStreamReader(in, "UTF-8"));
				try(PrintWriter out = super.openPathOrStdoutAsPrintWriter(this.outputFile)) {
				
					out.println("#chrom\tstart\tend"+(this.columns.isEmpty()?"":"\t")+String.join("\t",this.columns));

					String line;
					while((line=br.readLine())!=null) {
						if(line.startsWith("#") || StringUtils.isBlank(line)) continue;
						final String[] tokens = CharSplitter.TAB.split(line);
						if(tokens.length!=9) {
							throw new JvarkitException.TokenErrors(9,tokens);
						}
						
						final String ctg = contigNameConverter.apply(tokens[0]);
						if(StringUtils.isBlank(ctg)) continue;
						
						final Map<String,String> atts = attributes(tokens[8]);
						for(String pfx:fixed_cols_prefixes) {
							put(atts,pfx+"source",tokens[1]);
							put(atts,pfx+"feature",tokens[2]);
							put(atts,pfx+"score",tokens[5]);
							put(atts,pfx+"strand",tokens[6]);
							put(atts,pfx+"frame",tokens[7]);
						}
						
						
						boolean grep_is_ok=true;
						for(Grep g: greps) {
							if(!g.visit(atts)) grep_is_ok=false;
						}
						if(!grep_is_ok) continue;
						
						out.print(ctg);
						out.print('\t');
						out.print(Integer.parseInt(tokens[3])-1);
						out.print('\t');
						out.print(Integer.parseInt(tokens[4]));
						
						
						
						
						for(final String field:this.columns) {
							out.print('\t');
							out.print(atts.getOrDefault(field,"."));
							}
						out.println();
						}
					out.flush();
					}
				br.close();
				in.close();
				if(greps.stream().anyMatch(G->!G.isSuccess())) {
					LOG.error("end with error because one grep failed");
					return -1;
					}
				}
		
			return 0;
			}
		catch(final Throwable err) {
			LOG.error(err);
			return -1;
			}
		finally
			{
			}
		}
	
	public static void main(final String[] args)
		{
		new GtfToBed().instanceMainWithExit(args);
		}
	}
