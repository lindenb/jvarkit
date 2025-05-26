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
package com.github.lindenb.jvarkit.ucsc;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.OptionalInt;
import java.util.function.UnaryOperator;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import com.github.lindenb.jvarkit.io.FileHeader;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;
import com.github.lindenb.jvarkit.sql.parser.SchemaParser;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.readers.LineIterator;


public class UcscTranscriptCodec extends AsciiFeatureCodec<UcscTranscript> {
	private static final Logger LOG = Logger.of(UcscTranscriptCodec.class);

	
	static final String FILE_SUFFIX = ".txt.gz";
	private int contig_col =-1;
	private int name_col =-1;
	private int strand_col =-1;
	private int txStart_col =-1;
	private int txEnd_col =-1;
	private int cdsStart_col =-1;
	private int cdsEnd_col =-1;
	private int nExon_col =-1;
	private int exonsStart_col =-1;
	private int exonsEnd_col =-1;
	private int name2_col =-1;
	private int score_col =-1;
	private final boolean auto_detect_flag;
	private String sqlTableName="gene";
	private UnaryOperator<String> contigConverter = S->S;
	// file Header generating other columns
	private FileHeader fileHeader = null;
	
	public UcscTranscriptCodec() {
		super(UcscTranscript.class);
		auto_detect_flag =true;
		}
	
	public UcscTranscriptCodec(final SchemaParser.Table table ) {
		super(UcscTranscript.class);
		auto_detect_flag = false;
		sql(table);
		}
	public UcscTranscriptCodec(final Path p) {
		this(p.toString());
		}
	public UcscTranscriptCodec(final String uri) {
		super(UcscTranscript.class);
		auto_detect_flag = false;
		String sqluri;
		if(uri.endsWith(".sql")) {
			sqluri = uri;
		} else if(uri.endsWith(".txt") || uri.endsWith(".tsv")) {
			sqluri = uri.substring(0, uri.length() - 4) + ".sql";
		} else if(uri.endsWith(FILE_SUFFIX) || uri.endsWith(".tsv.gz")) {
			sqluri = uri.substring(0, uri.length() - FILE_SUFFIX.length()) + ".sql";
		} else
		{
			throw new IllegalArgumentException("uri must end with .sql or .txt or .tsv or .tsv.gz or "+FILE_SUFFIX);
			
		}
		try(InputStream in=IOUtils.openURIForReading(sqluri)) {
			final SchemaParser.Table table = SchemaParser.parseTable(in);
			sql(table);
		} catch(Throwable err) {
			throw new RuntimeIOException("cannot parse schema",err);
		}
	}
	
	/** get name of the SQL table associated to this codec */
	public String getName() {
		return sqlTableName;
		}
	
	private int columnIndex(final SchemaParser.Table table, String colName) {
		final SchemaParser.Column col = table.getColumnByName(colName).orElseThrow(()->new IllegalArgumentException("Cannot find column "+colName+" in table "+table.getName()));
		return col.getIndex();
		}
	
	private void sql(final SchemaParser.Table table ) {
		this.sqlTableName = table.getName();
		this.contig_col = columnIndex(table,"chrom");
		this.name_col = columnIndex(table,"name");
		this.strand_col = columnIndex(table,"strand");
		this.txStart_col = columnIndex(table,"txStart");
		this.txEnd_col = columnIndex(table,"txEnd");
		this.cdsStart_col = columnIndex(table,"cdsStart");
		this.cdsEnd_col = columnIndex(table,"cdsEnd");
		this.nExon_col = columnIndex(table,"exonCount");
		this.exonsStart_col = columnIndex(table,"exonStarts");
		this.exonsEnd_col = columnIndex(table,"exonEnds");
		if(table.getColumnByName("name2").isPresent()) {
			this.name2_col = columnIndex(table,"name2");
			}
		if(table.getColumnByName("name2").isPresent()) {
			this.score_col = columnIndex(table,"score");
			}
		this.fileHeader = new FileHeader(table.getColumns().stream().map(C->C.getName()).collect(Collectors.toList()));
		}
	
	public UcscTranscriptCodec setContigConverter(UnaryOperator<String> contigConverter) {
		this.contigConverter = contigConverter;
		return this;
		}
	@Override
	public UcscTranscript decode(final String s) {
		return decode(CharSplitter.TAB.split(s));
		}
	public UcscTranscript decode(final String[] s) {
		return decode(Arrays.asList(s));
		}
	public UcscTranscript decode(final List<String> tokens) {
		final UcscTranscriptImpl tr = new UcscTranscriptImpl();
		if(contig_col==-1) {
			if(auto_detect_flag) {
				/* first column may be the bin column 
				 * we use the strand to detect the format */
				final int binIdx=tokens.get(2).equals("+") || tokens.get(2).equals("-")?0:1;
				this.name_col = binIdx + 0;
				this.contig_col= binIdx + 1;
		        this.strand_col = binIdx + 2;
		        this.txStart_col = binIdx + 3;
		        this.txEnd_col = binIdx + 4;
		        this.cdsStart_col = binIdx + 5;
		        this.cdsEnd_col = binIdx + 6;
		        this.nExon_col = binIdx + 7;
		        this.exonsStart_col = binIdx + 8;
		        this.exonsEnd_col = binIdx + 9;
		        this.name2_col = binIdx + 11 ;
		        this.fileHeader= new FileHeader(IntStream.of(1,tokens.size()).mapToObj(i->"$"+i).collect(Collectors.toList()));
				}
			else
				{
				throw new IllegalStateException("undefined columns");
				}
			}
		tr.contig = contigConverter.apply(tokens.get(this.contig_col));
		if(StringUtils.isBlank(tr.contig)) return null;
		tr.name = tokens.get(this.name_col);
		tr.txStart = Integer.parseInt(tokens.get(this.txStart_col));
		tr.txEnd = Integer.parseInt(tokens.get(this.txEnd_col));
		String s = tokens.get(this.strand_col);
		if(s.equals("+")) {
			tr.strand = '+';
			}
		else if(s.equals("-")) {
			tr.strand = '-';
			}
		else
			{
			throw new IllegalArgumentException("bad strand in '"+s+"' for "+String.join(";", tokens));
			}
		tr.cdsStart = Integer.parseInt(tokens.get(this.cdsStart_col));
		tr.cdsEnd = Integer.parseInt(tokens.get(this.cdsEnd_col));
		final int nExons = Integer.parseInt(tokens.get(this.nExon_col));
		
		tr.exonStarts = new int[nExons];
		tr.exonEnds = new int[nExons];
		
		String[] ss = CharSplitter.COMMA.split(tokens.get(this.exonsStart_col));
		if(ss.length!=nExons) {
			throw new IllegalArgumentException("Expected "+nExons+" exonsStart but got "+ss.length+" for "+String.join(";", tokens));
			}
		for(int i=0;i< ss.length;i++) {
			tr.exonStarts[i] = Integer.parseInt(ss[i]);
			}
		
		ss = CharSplitter.COMMA.split(tokens.get(this.exonsEnd_col));
		if(ss.length!=nExons) {
			throw new IllegalArgumentException("Expected "+nExons+" exonsEnd but got "+ss.length+" for "+String.join(";", tokens));
			}
		for(int i=0;i< ss.length;i++) {
			tr.exonEnds[i] = Integer.parseInt(ss[i]);
			}
		if(this.name2_col>=0 && this.name2_col < tokens.size()) {
			tr.name2 = tokens.get(this.name2_col);
			}
		if(this.score_col>=0 && this.score_col < tokens.size()) {
			tr.score = OptionalInt.of(Integer.parseInt(tokens.get(this.score_col)));
			}
		
		// tr.metadata=this.fileHeader.toMap(tokens);
 		return tr;
		}
	
	private CloseableIterator<UcscTranscript> iterator(final String uri) throws IOException {
		return iterator(IOUtils.openURIForBufferedReading(uri));
		}

	
	private CloseableIterator<UcscTranscript> iterator(final Path p) throws IOException {
		return iterator(IOUtils.openPathForBufferedReading(p));
		}
	
	private CloseableIterator<UcscTranscript> iterator(final BufferedReader br) throws IOException {
		return new AbstractCloseableIterator<UcscTranscript>()
			{
			@Override
			protected UcscTranscript advance() {
				try {
					for(;;) {
						final String line = br.readLine();
						if(line==null) {
							close();
							return null;
							}
						if(StringUtils.isBlank(line) || line.startsWith("#")) continue;
						final UcscTranscript u = decode(line);
						if(u==null) continue;
						return u;
						}
					}
				catch(final IOException err) {
					throw new RuntimeIOException(err);
					}
				}
			@Override
			public void close() {
				try {br.close();}
				catch(Throwable err) {}
				}
			};
		}
	
	
	@Override
    public Object readActualHeader(final LineIterator lineIterator) {
        final List<String> header = new ArrayList<>();
		while (lineIterator.hasNext()) {
            final String nextLine = lineIterator.peek();
            if (nextLine.startsWith("#") || nextLine.isEmpty()) {
                // advance the iterator and consume the line (which is a no-op)
            	header.add(lineIterator.next());
            } else {
                return null; // break out when we've seen the end of the header
            }
        }

        return header;
    }

    @Override
    public boolean canDecode(final String path) {
        final String toDecode;
        if (IOUtil.hasBlockCompressedExtension(path)) {
            toDecode = path.substring(0, path.lastIndexOf("."));
        } else {
            toDecode = path;
        }
        return toDecode.endsWith(FILE_SUFFIX) || 
        		toDecode.endsWith(".txt") || 
        		toDecode.endsWith(".tsv") || 
        		toDecode.endsWith(".tsv.gz")
        		;
    }
    
    public static CloseableIterator<UcscTranscript> makeIterator(final InputStream in) throws IOException {
    	return makeIterator(in,null);
    	}
    
    public static CloseableIterator<UcscTranscript> makeIterator(final InputStream in,final String sqluri) throws IOException {
    	return makeIterator(new BufferedReader(new InputStreamReader(in)),sqluri);
    	}
    
    public static CloseableIterator<UcscTranscript> makeIterator(final BufferedReader br,final String sqluri) throws IOException {
    	final SchemaParser.Table table ;
    	if(StringUtils.isBlank(sqluri)) {
    		table = null;
    		}
    	else
    		{
    		table =  parseSchema(sqluri,true);
    		}
    	final UcscTranscriptCodec codec = table==null?new UcscTranscriptCodec():new UcscTranscriptCodec(table);
    	return codec.iterator(br);
    	}
    
    
    public static CloseableIterator<UcscTranscript> makeIterator(final String uri,final String sqluri) throws IOException {
    	final SchemaParser.Table table ;
    	if(StringUtils.isBlank(sqluri)) {
    		table = parseSchema(getSchemaUri(uri),false);
    		}
    	else
    		{
    		table =  parseSchema(sqluri,true);
    		}
    	return makeIterator(uri,table);
    	}
    	
    private static CloseableIterator<UcscTranscript> makeIterator(final String uri, final SchemaParser.Table schema) throws IOException {
    	final UcscTranscriptCodec codec = schema==null?new UcscTranscriptCodec():new UcscTranscriptCodec(schema);
    	return codec.iterator(uri);
    	}
    
    public static CloseableIterator<UcscTranscript> makeIterator(final String uri) throws IOException {
    	return makeIterator(uri,"");
    	}
    
    public static CloseableIterator<UcscTranscript> makeIterator(final Path path) throws IOException {
    	return makeIterator(path.toString(),"");
    	}

    
    private static String getSchemaUri(String uri) {
    	if(uri.endsWith(".gz")) uri =  uri.substring(0,uri.length()-3);
    	int dot= uri.lastIndexOf(".");
    	if(dot!=-1) uri=uri.substring(0,dot);
    	return uri+".sql";
    	}
    
    /** parse SQL schema at given uri */
    static SchemaParser.Table parseSchema(final String sqluri, boolean failOnError) throws IOException {
    	try(InputStream in=IOUtils.openURIForReading(sqluri)) {
			return SchemaParser.parseTable(in);
		} catch(Throwable err) {
			if(failOnError) throw new IOException("cannot parse sql schema from "+sqluri,err);
			LOG.warn("Cannot parse schema at "+sqluri);
			return null;
		}
    }
}
