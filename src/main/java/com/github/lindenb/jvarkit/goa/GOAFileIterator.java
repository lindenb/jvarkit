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
package com.github.lindenb.jvarkit.goa;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.iterator.AbstractCloseableIterator;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;

import htsjdk.samtools.util.RuntimeIOException;

/** GO Annotation File (GAF) files */
public class GOAFileIterator extends AbstractCloseableIterator<GOAFileIterator.GafRecord> {
	public static final String DEFAULT_GOA_URI ="http://cvsweb.geneontology.org/cgi-bin/cvsweb.cgi/go/gene-associations/gene_association.goa_human.gz?rev=HEAD";
	
	private final BufferedReader br;
	GOAFileIterator(final BufferedReader br) {
		this.br = br;
		}

	public interface GafRecord {
		/** eg UniProtKB */
		public String getDatabase();
		/** eg P12345 */
		public String getObjectId();
		/** eg PHO3 */
		public String getObjectSymbol();
		/** eg NOT */
		public Set<String> getQualifiers();
		/** eg GO:0003993 */
		public String getGoId();
		/** eg PMID:2676709 */
		public List<String> getDBReferences();
		/** eg IMP */
		public String getEvidenceCode();
		/** eg GO:0000346 */
		public Set<String> getWith();
		/** eg F */
		public String getAspect();
		/** eg Toll-like receptor 4  */
		public String getDBObjectName();
		/** eg   A0A024R5I4_HUMAN|SLCO2B1|hCG_27402 */
		public Set<String> getDBObjectSynonyms();
		/** A description of the type of gene product being annotated . Eg. protein */
		public String getDBObjectType();
		/** eg taxon:9606 */
		public Set<String> getTaxons();
		/** eg   */
		public String getDate();
		public String geAssigned();
		public String getAnnotation();
		public String getGeneProduct();
		}
	
	
	private class GafRecordImpl implements GafRecord {
		final List<String> tokens;
		
		GafRecordImpl(List<String> tokens) {
			this.tokens=tokens;
			}
		
		private String get(int i) {
			return this.tokens.get(i);
		}
		private String getMandatory(int i) {
			final String s= get(i);
			if(StringUtils.isBlank(s)) throw new IllegalStateException("missing field in column "+(i+1)+".");
			return s;
			}
		private Set<String> setOf(int i) {
			return new HashSet<>(listOf(i));
			}
		private List<String> listOf(int i) {
			final String s = get(i);
			if(StringUtils.isBlank(s)) return Collections.emptyList();
			return CharSplitter.PIPE.splitAsStringList(s);
		}
		/** eg UniProtKB */
		@Override public String getDatabase() { return getMandatory(0);}
		/** eg P12345 */
		@Override public String getObjectId() { return getMandatory(1);}
		/** eg PHO3 */
		@Override public String getObjectSymbol() { return get(2);}
		/** eg NOT */
		@Override public Set<String> getQualifiers() { return setOf(3);}
		/** eg GO:0003993 */
		@Override public String getGoId() { return getMandatory(4);}
		/** eg PMID:2676709 */
		@Override public List<String> getDBReferences() { return listOf(5);}
		/** eg IMP */
		@Override public String getEvidenceCode() { return get(6);}
		/** eg GO:0000346 */
		@Override public Set<String> getWith() { return setOf(7);}
		/** eg F */
		@Override public String getAspect() { return get(8);}
		/** eg Toll-like receptor 4  */
		@Override public String getDBObjectName() { return get(9);}
		/** eg   A0A024R5I4_HUMAN|SLCO2B1|hCG_27402 */
		@Override public Set<String> getDBObjectSynonyms() { return setOf(10);}
		/** A description of the type of gene product being annotated . Eg. protein */
		@Override public String getDBObjectType() { return get(11);}
		/** eg taxon:9606 */
		@Override public Set<String> getTaxons() { return setOf(12);}
		/** eg   */
		@Override public String getDate() { return get(13);}
		@Override public String geAssigned() { return get(14);}
		@Override public String getAnnotation() { return get(15);}
		@Override public String getGeneProduct() { return get(16);}
		
		@Override
		public String toString() {
			return String.join("\t", tokens);
			}
		}
	
	@Override
	protected GafRecord advance() {
		try {
			for(;;) {
				final String line = this.br.readLine();
				if(line==null) return null;
				if(line.startsWith("!")) continue;
				return new GafRecordImpl(CharSplitter.TAB.splitAsStringList(line));
				}
			} catch(final IOException err) {
				throw new RuntimeIOException(err);
				}
		}
	
	@Override
	public void close() {
		try {br.close();}
		catch(IOException err) {}
		}
	
	public static GOAFileIterator newInstance(String uri) throws IOException {
		return new GOAFileIterator(IOUtils.openURIForBufferedReading(uri));
		}
}
