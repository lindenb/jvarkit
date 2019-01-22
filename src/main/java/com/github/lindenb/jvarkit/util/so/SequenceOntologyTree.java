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
package com.github.lindenb.jvarkit.util.so;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.JvarkitException;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;



/**Sequence ontology tree */
public class SequenceOntologyTree
	implements Iterable<SequenceOntologyTree.Term>
	{
	private static final Logger LOG = Logger.build(SequenceOntologyTree.class).make();
	private static SequenceOntologyTree INSTANCE=null;
	private final Map<String,TermImpl> acn2term=new HashMap<>(3000);
	private final Map<String,TermImpl> label2term=new HashMap<>(3000);
	
	
	public interface Term
		{
		/** get unique uri for this term */
		public String getUri();
		/** get Term accession like SO:0000025 */
		public String getAcn();
		/** get Term label like "sugar_edge_base_pair" */
		public String getLabel();
		/** get the parent for this node */
		public Set<Term> getParents();
		/** get direct children of this node */
		public Set<Term> getChildren();
		/** get ALL (recursive) children of this node */
		public Set<Term> getAllDescendants();
		/** return true if term is children of parent */
		public boolean isChildrenOf(final Term t);
		}
	private class TermImpl implements Term
		{
		final String accession;
		final int _hash;
		String label;
		final Set<Term> parents=new HashSet<>();
		final Set<Term> children=new HashSet<>();
		
		TermImpl(final String accession,final String label) {
			this.accession = accession;
			this._hash = accession.hashCode();
			this.label = label;/* may be null */
		}
		/** get URL "http://purl.obolibrary.org/obo/..."  */
		@Override
		public String getUri() {
			return "http://purl.obolibrary.org/obo/"+getAcn().replace(':', '_');
			}
		
		public String getAcn() {
			return accession;
			}
		@Override
		public String getLabel() {
			return label;
			}
		
		
		/** returns only the direct children of this node */
		@Override
		public Set<Term> getChildren()
			{
			return Collections.unmodifiableSet(this.children);
			}
		
		/** recursive operation on getChildren, including self */
		@Override
		public Set<Term> getAllDescendants()
			{
			final Set<Term> set=new HashSet<Term>();
			_getAllDescendants(this,set);
			return set;
			}
		/** return true if term is children of parent */
		public boolean isChildrenOf(final Term t) {
			return _isChildrenOf(this,t);
			}	

		
		@Override
		public Set<Term> getParents()
			{
			return Collections.unmodifiableSet(this.parents);
			}
		
		@Override
		public int hashCode()
			{
			return this._hash;
			}
		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			final TermImpl other = (TermImpl) obj;
			return this.accession.equals(other.accession);
			}
		@Override
		public String toString() {
			return accession;
			}

		}
	

	
	private static void _getAllDescendants(final Term term,final Set<Term> set)
		{
		set.add(term);
		for(final Term c:term.getChildren())
			{
			_getAllDescendants(c,set);
			}
		}

	private static boolean _isChildrenOf(final Term term,final Term parent)
		{
		if(parent.equals(term)) return true;
		for(final Term c:parent.getChildren())
			{
			if(_isChildrenOf(term,c)) return true;
			}
		return false;
		}
	
	/*
	private TermImpl createTerm(final String acn,final String label)
		{
		final TermImpl t=new TermImpl(acn,label);
		this.acn2term.put(t.accession, t);
		this.label2term.put(t.getLabel(), t);
		return t;
		}*/
	
	private TermImpl newTerm(final String acn,final String label,final String...parentAcns)
		{
		TermImpl term = this.acn2term.get(acn);
		if(term ==null) {
			term =new TermImpl(acn, label);
			this.acn2term.put(acn, term);
			this.label2term.put(label, term);
			}
		else if(term.label==null)
			{
			term.label = label;
			this.label2term.put(label, term);
			}
		else {
			throw new IllegalArgumentException("term "+acn+" label ");
			}
		for(final String parentAcn: parentAcns) {
			TermImpl parentTerm = this.acn2term.get(parentAcn);
			if(parentTerm==null) {
				parentTerm =new TermImpl(parentAcn, null);
				this.acn2term.put(parentAcn, parentTerm);
				}
			parentTerm.children.add(term);
			term.parents.add(parentTerm);
			}
		return term;
		}
	
	public Term getTermByAcn(final String s)
		{
		return this.acn2term.get(s);
		}
	
	private static String normalizeName(final String s)
		{
		return s.toLowerCase().
				replaceAll("'", "prime").
				replaceAll(" ", "_");
		}
	
	/** loop over terms and find a term.label==user.label */ 
	public Term getTermByLabel(final String s)
		{
		final Term t= this.label2term.get(s);
		if(t!=null) return t;
		return  this.label2term.get(normalizeName(s));
		}
	
	private SequenceOntologyTree()
		{
		
		}
	/** get all terms available in this tree */
	public Collection<? extends Term> getTerms()
		{
		return Collections.unmodifiableCollection(this.acn2term.values());
		}
	
	@SuppressWarnings("unchecked")
	@Override
	public Iterator<Term> iterator() {
		return ( Iterator<Term>)getTerms().iterator();
		}
	
	public String getLastUpdated()
		{
		return "2017-03-31";
		}
	
	public static  SequenceOntologyTree fromUri(final String uri) throws IOException {
		return new OwlLoader().parse(uri);
		}

	
	
	public static  SequenceOntologyTree getInstance()
		{
		if(INSTANCE==null)
			{
			synchronized (SequenceOntologyTree.class)
				{
				if(INSTANCE==null)
					{
					INSTANCE= createNewInstance();
					}
				}
			}
		return INSTANCE;
		}
	
	/** try first to read from /META-INF/so/so-simple.owl using createFromInputStream if it fails calls createDefault */
	public static SequenceOntologyTree createNewInstance() {
		try(final InputStream in =SequenceOntologyTree.class.getResourceAsStream("/META-INF/so/so-simple.owl");)
			{
			if(in!=null) return SequenceOntologyTree.createFromInputStream(in);
			}
		catch(final Exception error) {
			// ignore
			}
		return SequenceOntologyTree.createDefault();	
		}

	public static SequenceOntologyTree createFromInputStream(final InputStream in) throws IOException
		{
		return new OwlLoader().parse(in);
		}
	public static SequenceOntologyTree createDefault()
		{
		final SequenceOntologyTree tree =new SequenceOntologyTree();
		//curl -L  "https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/master/so-simple.owl?raw=true" | xsltproc so2java.xsl -
		tree.newTerm("SO:0000001","region","SO:0000110");
		tree.newTerm("SO:0000002","sequence_secondary_structure","SO:0001411");
		tree.newTerm("SO:0000003","G_quartet","SO:0000002");
		tree.newTerm("SO:0000004","interior_coding_exon","SO:0000195");
		tree.newTerm("SO:0000005","satellite_DNA","SO:0000705");
		tree.newTerm("SO:0000006","PCR_product","SO:0000695");
		tree.newTerm("SO:0000007","read_pair","SO:0000150");
		tree.newTerm("SO:0000010","protein_coding","SO:0000401");
		tree.newTerm("SO:0000011","non_protein_coding","SO:0000401");
		tree.newTerm("SO:0000012","scRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000013","scRNA","SO:0000655");
		tree.newTerm("SO:0000014","INR_motif","SO:0001660");
		tree.newTerm("SO:0000015","DPE_motif","SO:0001660");
		tree.newTerm("SO:0000016","BREu_motif","SO:0001660");
		tree.newTerm("SO:0000017","PSE_motif","SO:0000713");
		tree.newTerm("SO:0000018","linkage_group","SO:0001411");
		tree.newTerm("SO:0000020","RNA_internal_loop","SO:0000715");
		tree.newTerm("SO:0000021","asymmetric_RNA_internal_loop","SO:0000020");
		tree.newTerm("SO:0000022","A_minor_RNA_motif","SO:0000715");
		tree.newTerm("SO:0000023","K_turn_RNA_motif","SO:0000021");
		tree.newTerm("SO:0000024","sarcin_like_RNA_motif","SO:0000021");
		tree.newTerm("SO:0000025","symmetric_RNA_internal_loop","SO:0000020");
		tree.newTerm("SO:0000026","RNA_junction_loop","SO:0000715");
		tree.newTerm("SO:0000027","RNA_hook_turn","SO:0000026");
		tree.newTerm("SO:0000028","base_pair","SO:0000002");
		tree.newTerm("SO:0000029","WC_base_pair","SO:0000028");
		tree.newTerm("SO:0000030","sugar_edge_base_pair","SO:0000028");
		tree.newTerm("SO:0000031","aptamer","SO:0000696");
		tree.newTerm("SO:0000032","DNA_aptamer","SO:0000031");
		tree.newTerm("SO:0000033","RNA_aptamer","SO:0000031");
		tree.newTerm("SO:0000034","morpholino_oligo","SO:0001247");
		tree.newTerm("SO:0000035","riboswitch","SO:0000836");
		tree.newTerm("SO:0000036","matrix_attachment_site","SO:0000626");
		tree.newTerm("SO:0000037","locus_control_region","SO:0000727");
		tree.newTerm("SO:0000039","match_part","SO:0001410");
		tree.newTerm("SO:0000040","genomic_clone","SO:0000151");
		tree.newTerm("SO:0000043","processed_pseudogene","SO:0000336");
		tree.newTerm("SO:0000044","pseudogene_by_unequal_crossing_over","SO:0001760");
		tree.newTerm("SO:0000051","probe","SO:0000696");
		tree.newTerm("SO:0000054","aneuploid","SO:1000182");
		tree.newTerm("SO:0000055","hyperploid","SO:0000054");
		tree.newTerm("SO:0000056","hypoploid","SO:0000054");
		tree.newTerm("SO:0000057","operator","SO:0000752");
		tree.newTerm("SO:0000059","nuclease_binding_site","SO:0001654");
		tree.newTerm("SO:0000060","compound_chromosome_arm","SO:1000042");
		tree.newTerm("SO:0000061","restriction_enzyme_binding_site","SO:0000059");
		tree.newTerm("SO:0000062","deficient_intrachromosomal_transposition","SO:1000029","SO:1000041");
		tree.newTerm("SO:0000063","deficient_interchromosomal_transposition","SO:1000155");
		tree.newTerm("SO:0000065","free_chromosome_arm","SO:1000183");
		tree.newTerm("SO:0000067","gene_to_gene_feature","SO:0000401");
		tree.newTerm("SO:0000068","overlapping","SO:0000067");
		tree.newTerm("SO:0000069","inside_intron","SO:0000068");
		tree.newTerm("SO:0000070","inside_intron_antiparallel","SO:0000069");
		tree.newTerm("SO:0000071","inside_intron_parallel","SO:0000069");
		tree.newTerm("SO:0000073","five_prime_three_prime_overlap","SO:0000068");
		tree.newTerm("SO:0000074","five_prime_five_prime_overlap","SO:0000068");
		tree.newTerm("SO:0000075","three_prime_three_prime_overlap","SO:0000068");
		tree.newTerm("SO:0000076","three_prime_five_prime_overlap","SO:0000068");
		tree.newTerm("SO:0000077","antisense","SO:0000068");
		tree.newTerm("SO:0000078","polycistronic_transcript","SO:0000673");
		tree.newTerm("SO:0000079","dicistronic_transcript","SO:0000078");
		tree.newTerm("SO:0000080","operon_member","SO:0000081");
		tree.newTerm("SO:0000081","gene_array_member","SO:0000401");
		tree.newTerm("SO:0000083","macronuclear_sequence","SO:0000736");
		tree.newTerm("SO:0000084","micronuclear_sequence","SO:0000736");
		tree.newTerm("SO:0000087","nuclear_gene","SO:0000704");
		tree.newTerm("SO:0000088","mt_gene","SO:0000704");
		tree.newTerm("SO:0000089","kinetoplast_gene","SO:0000088");
		tree.newTerm("SO:0000090","plastid_gene","SO:0000704");
		tree.newTerm("SO:0000091","apicoplast_gene","SO:0000090");
		tree.newTerm("SO:0000092","ct_gene","SO:0000090");
		tree.newTerm("SO:0000093","chromoplast_gene","SO:0000090");
		tree.newTerm("SO:0000094","cyanelle_gene","SO:0000090");
		tree.newTerm("SO:0000095","leucoplast_gene","SO:0000090");
		tree.newTerm("SO:0000096","proplastid_gene","SO:0000090");
		tree.newTerm("SO:0000097","nucleomorph_gene","SO:0000704");
		tree.newTerm("SO:0000098","plasmid_gene","SO:0000704");
		tree.newTerm("SO:0000099","proviral_gene","SO:0000704");
		tree.newTerm("SO:0000100","endogenous_retroviral_gene","SO:0000099");
		tree.newTerm("SO:0000101","transposable_element","SO:0001039");
		tree.newTerm("SO:0000102","expressed_sequence_match","SO:0000347");
		tree.newTerm("SO:0000103","clone_insert_end","SO:0000699");
		tree.newTerm("SO:0000104","polypeptide","SO:0001411");
		tree.newTerm("SO:0000105","chromosome_arm","SO:0000830");
		tree.newTerm("SO:0000107","sequencing_primer","SO:0000112");
		tree.newTerm("SO:0000108","mRNA_with_frameshift","SO:0000234");
		tree.newTerm("SO:0000110","sequence_feature");
		tree.newTerm("SO:0000111","transposable_element_gene","SO:0000704");
		tree.newTerm("SO:0000112","primer","SO:0000441");
		tree.newTerm("SO:0000113","proviral_region","SO:0001039");
		tree.newTerm("SO:0000114","methylated_cytosine","SO:0000306","SO:0001963");
		tree.newTerm("SO:0000116","edited","SO:0000237");
		tree.newTerm("SO:0000118","transcript_with_translational_frameshift","SO:0000673");
		tree.newTerm("SO:0000119","regulated","SO:0000401");
		tree.newTerm("SO:0000120","protein_coding_primary_transcript","SO:0000185");
		tree.newTerm("SO:0000121","forward_primer","SO:0000112");
		tree.newTerm("SO:0000122","RNA_sequence_secondary_structure","SO:0000002");
		tree.newTerm("SO:0000123","transcriptionally_regulated","SO:0000119");
		tree.newTerm("SO:0000124","transcriptionally_constitutive","SO:0000123");
		tree.newTerm("SO:0000125","transcriptionally_induced","SO:0000123");
		tree.newTerm("SO:0000126","transcriptionally_repressed","SO:0000123");
		tree.newTerm("SO:0000127","silenced_gene","SO:0000704");
		tree.newTerm("SO:0000128","gene_silenced_by_DNA_modification","SO:0000127");
		tree.newTerm("SO:0000129","gene_silenced_by_DNA_methylation","SO:0000128");
		tree.newTerm("SO:0000130","post_translationally_regulated","SO:0000119");
		tree.newTerm("SO:0000131","translationally_regulated","SO:0000119");
		tree.newTerm("SO:0000132","reverse_primer","SO:0000112");
		tree.newTerm("SO:0000133","epigenetically_modified","SO:0000401");
		tree.newTerm("SO:0000134","genomically_imprinted","SO:0000119","SO:0000133");
		tree.newTerm("SO:0000135","maternally_imprinted","SO:0000134");
		tree.newTerm("SO:0000136","paternally_imprinted","SO:0000134");
		tree.newTerm("SO:0000137","allelically_excluded","SO:0000133");
		tree.newTerm("SO:0000138","gene_rearranged_at_DNA_level","SO:0000898");
		tree.newTerm("SO:0000139","ribosome_entry_site","SO:0000836");
		tree.newTerm("SO:0000140","attenuator","SO:0001680");
		tree.newTerm("SO:0000141","terminator","SO:0001679");
		tree.newTerm("SO:0000142","DNA_sequence_secondary_structure","SO:0000002");
		tree.newTerm("SO:0000143","assembly_component","SO:0001410");
		tree.newTerm("SO:0000145","recoded_codon","SO:0000360");
		tree.newTerm("SO:0000146","capped","SO:0000237");
		tree.newTerm("SO:0000147","exon","SO:0000833");
		tree.newTerm("SO:0000148","supercontig","SO:0001876");
		tree.newTerm("SO:0000149","contig","SO:0000143","SO:0000353");
		tree.newTerm("SO:0000150","read","SO:0000143");
		tree.newTerm("SO:0000151","clone","SO:0000695");
		tree.newTerm("SO:0000152","YAC","SO:0000440");
		tree.newTerm("SO:0000153","BAC","SO:0000440");
		tree.newTerm("SO:0000154","PAC","SO:0000440");
		tree.newTerm("SO:0000155","plasmid","SO:0001235");
		tree.newTerm("SO:0000156","cosmid","SO:0000440");
		tree.newTerm("SO:0000157","phagemid","SO:0000440");
		tree.newTerm("SO:0000158","fosmid","SO:0000440");
		tree.newTerm("SO:0000159","deletion","SO:0001059","SO:0001411");
		tree.newTerm("SO:0000161","methylated_adenine","SO:0000306","SO:0001962");
		tree.newTerm("SO:0000162","splice_site","SO:0000835");
		tree.newTerm("SO:0000163","five_prime_cis_splice_site","SO:0001419");
		tree.newTerm("SO:0000164","three_prime_cis_splice_site","SO:0001419");
		tree.newTerm("SO:0000165","enhancer","SO:0000727");
		tree.newTerm("SO:0000166","enhancer_bound_by_factor","SO:0000165");
		tree.newTerm("SO:0000167","promoter","SO:0001055");
		tree.newTerm("SO:0000169","RNApol_I_promoter","SO:0001203");
		tree.newTerm("SO:0000170","RNApol_II_promoter","SO:0001203");
		tree.newTerm("SO:0000171","RNApol_III_promoter","SO:0001203");
		tree.newTerm("SO:0000172","CAAT_signal","SO:0000713");
		tree.newTerm("SO:0000173","GC_rich_promoter_region","SO:0001659");
		tree.newTerm("SO:0000174","TATA_box","SO:0001660");
		tree.newTerm("SO:0000175","minus_10_signal","SO:0000713");
		tree.newTerm("SO:0000176","minus_35_signal","SO:0000713");
		tree.newTerm("SO:0000177","cross_genome_match","SO:0000347");
		tree.newTerm("SO:0000178","operon","SO:0005855");
		tree.newTerm("SO:0000179","clone_insert_start","SO:0000699");
		tree.newTerm("SO:0000180","retrotransposon","SO:0000101");
		tree.newTerm("SO:0000181","translated_nucleotide_match","SO:0000347");
		tree.newTerm("SO:0000182","DNA_transposon","SO:0000101");
		tree.newTerm("SO:0000183","non_transcribed_region","SO:0000842");
		tree.newTerm("SO:0000184","U2_intron","SO:0000662");
		tree.newTerm("SO:0000185","primary_transcript","SO:0000673");
		tree.newTerm("SO:0000186","LTR_retrotransposon","SO:0000180");
		tree.newTerm("SO:0000188","intron","SO:0000835");
		tree.newTerm("SO:0000189","non_LTR_retrotransposon","SO:0000180");
		tree.newTerm("SO:0000190","five_prime_intron","SO:0000188");
		tree.newTerm("SO:0000191","interior_intron","SO:0000188");
		tree.newTerm("SO:0000192","three_prime_intron","SO:0000188");
		tree.newTerm("SO:0000193","RFLP_fragment","SO:0000412");
		tree.newTerm("SO:0000194","LINE_element","SO:0000189");
		tree.newTerm("SO:0000195","coding_exon","SO:0000147");
		tree.newTerm("SO:0000196","five_prime_coding_exon_coding_region","SO:0001215");
		tree.newTerm("SO:0000197","three_prime_coding_exon_coding_region","SO:0001215");
		tree.newTerm("SO:0000198","noncoding_exon","SO:0000147");
		tree.newTerm("SO:0000199","translocation","SO:0001785");
		tree.newTerm("SO:0000200","five_prime_coding_exon","SO:0000195");
		tree.newTerm("SO:0000201","interior_exon","SO:0000147");
		tree.newTerm("SO:0000202","three_prime_coding_exon","SO:0000195");
		tree.newTerm("SO:0000203","UTR","SO:0000836");
		tree.newTerm("SO:0000204","five_prime_UTR","SO:0000203");
		tree.newTerm("SO:0000205","three_prime_UTR","SO:0000203");
		tree.newTerm("SO:0000206","SINE_element","SO:0000189");
		tree.newTerm("SO:0000207","simple_sequence_length_variation","SO:0000248");
		tree.newTerm("SO:0000208","terminal_inverted_repeat_element","SO:0000182");
		tree.newTerm("SO:0000209","rRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000210","tRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000211","alanine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000212","arginine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000213","asparagine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000214","aspartic_acid_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000215","cysteine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000216","glutamic_acid_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000217","glutamine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000218","glycine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000219","histidine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000220","isoleucine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000221","leucine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000222","lysine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000223","methionine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000224","phenylalanine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000225","proline_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000226","serine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000227","threonine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000228","tryptophan_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000229","tyrosine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000230","valine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0000231","snRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000232","snoRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000233","mature_transcript","SO:0000673");
		tree.newTerm("SO:0000234","mRNA","SO:0000233");
		tree.newTerm("SO:0000235","TF_binding_site","SO:0001654","SO:0001679");
		tree.newTerm("SO:0000236","ORF","SO:0000717");
		tree.newTerm("SO:0000237","transcript_attribute","SO:0000733");
		tree.newTerm("SO:0000238","foldback_element","SO:0000182");
		tree.newTerm("SO:0000239","flanking_region","SO:0001412");
		tree.newTerm("SO:0000240","chromosome_variation","SO:0001507");
		tree.newTerm("SO:0000241","internal_UTR","SO:0000203");
		tree.newTerm("SO:0000242","untranslated_region_polycistronic_mRNA","SO:0000203");
		tree.newTerm("SO:0000243","internal_ribosome_entry_site","SO:0000139");
		tree.newTerm("SO:0000246","polyadenylated","SO:0000863");
		tree.newTerm("SO:0000248","sequence_length_variation","SO:0001059");
		tree.newTerm("SO:0000250","modified_RNA_base_feature","SO:0001236");
		tree.newTerm("SO:0000252","rRNA","SO:0000655");
		tree.newTerm("SO:0000253","tRNA","SO:0000655");
		tree.newTerm("SO:0000254","alanyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000255","rRNA_small_subunit_primary_transcript","SO:0000209");
		tree.newTerm("SO:0000256","asparaginyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000257","aspartyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000258","cysteinyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000259","glutaminyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000260","glutamyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000261","glycyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000262","histidyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000263","isoleucyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000264","leucyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000265","lysyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000266","methionyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000267","phenylalanyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000268","prolyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000269","seryl_tRNA","SO:0000253");
		tree.newTerm("SO:0000270","threonyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000271","tryptophanyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000272","tyrosyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000273","valyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000274","snRNA","SO:0000655");
		tree.newTerm("SO:0000275","snoRNA","SO:0000655");
		tree.newTerm("SO:0000276","miRNA","SO:0000370");
		tree.newTerm("SO:0000277","bound_by_factor","SO:0000733");
		tree.newTerm("SO:0000278","transcript_bound_by_nucleic_acid","SO:0000673");
		tree.newTerm("SO:0000279","transcript_bound_by_protein","SO:0000673");
		tree.newTerm("SO:0000280","engineered_gene","SO:0000704","SO:0000804");
		tree.newTerm("SO:0000281","engineered_foreign_gene","SO:0000280","SO:0000285","SO:0000805");
		tree.newTerm("SO:0000282","mRNA_with_minus_1_frameshift","SO:0000108");
		tree.newTerm("SO:0000283","engineered_foreign_transposable_element_gene","SO:0000111","SO:0000281");
		tree.newTerm("SO:0000285","foreign_gene","SO:0000704");
		tree.newTerm("SO:0000286","long_terminal_repeat","SO:0000657");
		tree.newTerm("SO:0000287","fusion_gene","SO:0000704");
		tree.newTerm("SO:0000288","engineered_fusion_gene","SO:0000280","SO:0000287");
		tree.newTerm("SO:0000289","microsatellite","SO:0000005");
		tree.newTerm("SO:0000290","dinucleotide_repeat_microsatellite_feature","SO:0000289");
		tree.newTerm("SO:0000291","trinucleotide_repeat_microsatellite_feature","SO:0000289");
		tree.newTerm("SO:0000293","engineered_foreign_repetitive_element","SO:0000657","SO:0000805");
		tree.newTerm("SO:0000294","inverted_repeat","SO:0000657");
		tree.newTerm("SO:0000295","U12_intron","SO:0000662");
		tree.newTerm("SO:0000296","origin_of_replication","SO:0001411");
		tree.newTerm("SO:0000297","D_loop","SO:0000296");
		tree.newTerm("SO:0000298","recombination_feature","SO:0001411");
		tree.newTerm("SO:0000299","specific_recombination_site","SO:0000669");
		tree.newTerm("SO:0000300","recombination_feature_of_rearranged_gene","SO:0000299");
		tree.newTerm("SO:0000301","vertebrate_immune_system_gene_recombination_feature","SO:0000300");
		tree.newTerm("SO:0000302","J_gene_recombination_feature","SO:0000939");
		tree.newTerm("SO:0000303","clip","SO:0000835");
		tree.newTerm("SO:0000305","modified_DNA_base","SO:0001236","SO:0001720");
		tree.newTerm("SO:0000306","methylated_DNA_base_feature","SO:0000305");
		tree.newTerm("SO:0000307","CpG_island","SO:0001411");
		tree.newTerm("SO:0000312","experimentally_determined","SO:0000789");
		tree.newTerm("SO:0000313","stem_loop","SO:0000122");
		tree.newTerm("SO:0000314","direct_repeat","SO:0000657");
		tree.newTerm("SO:0000315","TSS","SO:0000835");
		tree.newTerm("SO:0000316","CDS","SO:0000836");
		tree.newTerm("SO:0000317","cDNA_clone","SO:0000151");
		tree.newTerm("SO:0000318","start_codon","SO:0000360");
		tree.newTerm("SO:0000319","stop_codon","SO:0000360");
		tree.newTerm("SO:0000320","intronic_splice_enhancer","SO:0000344");
		tree.newTerm("SO:0000321","mRNA_with_plus_1_frameshift","SO:0000108");
		tree.newTerm("SO:0000322","nuclease_hypersensitive_site","SO:0000684");
		tree.newTerm("SO:0000323","coding_start","SO:0000851");
		tree.newTerm("SO:0000324","tag","SO:0000696");
		tree.newTerm("SO:0000325","rRNA_large_subunit_primary_transcript","SO:0000209");
		tree.newTerm("SO:0000326","SAGE_tag","SO:0000324");
		tree.newTerm("SO:0000327","coding_end","SO:0000851");
		tree.newTerm("SO:0000328","microarray_oligo","SO:0000051");
		tree.newTerm("SO:0000329","mRNA_with_plus_2_frameshift","SO:0000108");
		tree.newTerm("SO:0000330","conserved_region","SO:0001410");
		tree.newTerm("SO:0000331","STS","SO:0000324");
		tree.newTerm("SO:0000332","coding_conserved_region","SO:0000330");
		tree.newTerm("SO:0000333","exon_junction","SO:0000699");
		tree.newTerm("SO:0000334","nc_conserved_region","SO:0000330");
		tree.newTerm("SO:0000335","mRNA_with_minus_2_frameshift","SO:0000108");
		tree.newTerm("SO:0000336","pseudogene","SO:0001411");
		tree.newTerm("SO:0000337","RNAi_reagent","SO:0000442");
		tree.newTerm("SO:0000338","MITE","SO:0000208");
		tree.newTerm("SO:0000339","recombination_hotspot","SO:0000298");
		tree.newTerm("SO:0000340","chromosome","SO:0001235");
		tree.newTerm("SO:0000341","chromosome_band","SO:0000830");
		tree.newTerm("SO:0000342","site_specific_recombination_target_region","SO:0000299");
		tree.newTerm("SO:0000343","match","SO:0001410");
		tree.newTerm("SO:0000344","splice_enhancer","SO:0000165","SO:0001056");
		tree.newTerm("SO:0000345","EST","SO:0000324");
		tree.newTerm("SO:0000346","loxP_site","SO:0000947");
		tree.newTerm("SO:0000347","nucleotide_match","SO:0000343");
		tree.newTerm("SO:0000348","nucleic_acid","SO:0000443");
		tree.newTerm("SO:0000349","protein_match","SO:0000343");
		tree.newTerm("SO:0000350","FRT_site","SO:0000948");
		tree.newTerm("SO:0000351","synthetic_sequence","SO:0000443");
		tree.newTerm("SO:0000352","DNA","SO:0000348");
		tree.newTerm("SO:0000353","sequence_assembly","SO:0001248");
		tree.newTerm("SO:0000354","group_1_intron_homing_endonuclease_target_region","SO:0000684");
		tree.newTerm("SO:0000355","haplotype_block","SO:0000298");
		tree.newTerm("SO:0000356","RNA","SO:0000348");
		tree.newTerm("SO:0000357","flanked","SO:0000733");
		tree.newTerm("SO:0000359","floxed","SO:0000357");
		tree.newTerm("SO:0000360","codon","SO:0000851");
		tree.newTerm("SO:0000361","FRT_flanked","SO:0000357");
		tree.newTerm("SO:0000362","invalidated_by_chimeric_cDNA","SO:0000790");
		tree.newTerm("SO:0000363","floxed_gene","SO:0000902");
		tree.newTerm("SO:0000364","transposable_element_flanking_region","SO:0000239");
		tree.newTerm("SO:0000365","integron","SO:0001039");
		tree.newTerm("SO:0000366","insertion_site","SO:0000699");
		tree.newTerm("SO:0000367","attI_site","SO:0000946");
		tree.newTerm("SO:0000368","transposable_element_insertion_site","SO:0000366");
		tree.newTerm("SO:0000370","small_regulatory_ncRNA","SO:0000655");
		tree.newTerm("SO:0000371","conjugative_transposon","SO:0000182");
		tree.newTerm("SO:0000372","enzymatic_RNA","SO:0000673");
		tree.newTerm("SO:0000373","recombinationally_inverted_gene","SO:0000456");
		tree.newTerm("SO:0000374","ribozyme","SO:0000372");
		tree.newTerm("SO:0000375","rRNA_5_8S","SO:0000651");
		tree.newTerm("SO:0000376","RNA_6S","SO:0000370");
		tree.newTerm("SO:0000377","CsrB_RsmB_RNA","SO:0000370");
		tree.newTerm("SO:0000378","DsrA_RNA","SO:0000370");
		tree.newTerm("SO:0000379","GcvB_RNA","SO:0000378");
		tree.newTerm("SO:0000380","hammerhead_ribozyme","SO:0000715");
		tree.newTerm("SO:0000381","group_IIA_intron","SO:0000603");
		tree.newTerm("SO:0000382","group_IIB_intron","SO:0000603");
		tree.newTerm("SO:0000383","MicF_RNA","SO:0000644");
		tree.newTerm("SO:0000384","OxyS_RNA","SO:0000370");
		tree.newTerm("SO:0000385","RNase_MRP_RNA","SO:0000655");
		tree.newTerm("SO:0000386","RNase_P_RNA","SO:0000655");
		tree.newTerm("SO:0000387","RprA_RNA","SO:0000370");
		tree.newTerm("SO:0000388","RRE_RNA","SO:0000370");
		tree.newTerm("SO:0000389","spot_42_RNA","SO:0000370");
		tree.newTerm("SO:0000390","telomerase_RNA","SO:0000655");
		tree.newTerm("SO:0000391","U1_snRNA","SO:0000274");
		tree.newTerm("SO:0000392","U2_snRNA","SO:0000274");
		tree.newTerm("SO:0000393","U4_snRNA","SO:0000274");
		tree.newTerm("SO:0000394","U4atac_snRNA","SO:0000274");
		tree.newTerm("SO:0000395","U5_snRNA","SO:0000274");
		tree.newTerm("SO:0000396","U6_snRNA","SO:0000274");
		tree.newTerm("SO:0000397","U6atac_snRNA","SO:0000274");
		tree.newTerm("SO:0000398","U11_snRNA","SO:0000274");
		tree.newTerm("SO:0000399","U12_snRNA","SO:0000274");
		tree.newTerm("SO:0000400","sequence_attribute");
		tree.newTerm("SO:0000401","gene_attribute","SO:0000733");
		tree.newTerm("SO:0000403","U14_snoRNA","SO:0000593");
		tree.newTerm("SO:0000404","vault_RNA","SO:0000655");
		tree.newTerm("SO:0000405","Y_RNA","SO:0000655");
		tree.newTerm("SO:0000406","twintron","SO:0000188");
		tree.newTerm("SO:0000407","rRNA_18S","SO:0000650");
		tree.newTerm("SO:0000409","binding_site","SO:0001411");
		tree.newTerm("SO:0000410","protein_binding_site","SO:0000409");
		tree.newTerm("SO:0000411","rescue_region","SO:0000695");
		tree.newTerm("SO:0000412","restriction_fragment","SO:0000143");
		tree.newTerm("SO:0000413","sequence_difference","SO:0000700");
		tree.newTerm("SO:0000414","invalidated_by_genomic_contamination","SO:0000790");
		tree.newTerm("SO:0000415","invalidated_by_genomic_polyA_primed_cDNA","SO:0000790");
		tree.newTerm("SO:0000416","invalidated_by_partial_processing","SO:0000790");
		tree.newTerm("SO:0000417","polypeptide_domain","SO:0001070","SO:0100021");
		tree.newTerm("SO:0000418","signal_peptide","SO:0001527");
		tree.newTerm("SO:0000419","mature_protein_region","SO:0000839");
		tree.newTerm("SO:0000420","five_prime_terminal_inverted_repeat","SO:0000481");
		tree.newTerm("SO:0000421","three_prime_terminal_inverted_repeat","SO:0000481");
		tree.newTerm("SO:0000422","U5_LTR_region","SO:0000848");
		tree.newTerm("SO:0000423","R_LTR_region","SO:0000848");
		tree.newTerm("SO:0000424","U3_LTR_region","SO:0000848");
		tree.newTerm("SO:0000425","five_prime_LTR","SO:0000286");
		tree.newTerm("SO:0000426","three_prime_LTR","SO:0000286");
		tree.newTerm("SO:0000427","R_five_prime_LTR_region","SO:0000423","SO:0000850");
		tree.newTerm("SO:0000428","U5_five_prime_LTR_region","SO:0000422","SO:0000850");
		tree.newTerm("SO:0000429","U3_five_prime_LTR_region","SO:0000424","SO:0000850");
		tree.newTerm("SO:0000430","R_three_prime_LTR_region","SO:0000849");
		tree.newTerm("SO:0000431","U3_three_prime_LTR_region","SO:0000849");
		tree.newTerm("SO:0000432","U5_three_prime_LTR_region","SO:0000849");
		tree.newTerm("SO:0000433","non_LTR_retrotransposon_polymeric_tract","SO:0000657","SO:0000840");
		tree.newTerm("SO:0000434","target_site_duplication","SO:0000314");
		tree.newTerm("SO:0000435","RR_tract","SO:0000330");
		tree.newTerm("SO:0000436","ARS","SO:0000296");
		tree.newTerm("SO:0000439","inverted_ring_chromosome","SO:1000030","SO:1000045");
		tree.newTerm("SO:0000440","vector_replicon","SO:0001235");
		tree.newTerm("SO:0000441","ss_oligo","SO:0000696");
		tree.newTerm("SO:0000442","ds_oligo","SO:0000696");
		tree.newTerm("SO:0000443","polymer_attribute","SO:0000400");
		tree.newTerm("SO:0000444","three_prime_noncoding_exon","SO:0000198");
		tree.newTerm("SO:0000445","five_prime_noncoding_exon","SO:0000198");
		tree.newTerm("SO:0000446","UTR_intron","SO:0000188");
		tree.newTerm("SO:0000447","five_prime_UTR_intron","SO:0000446");
		tree.newTerm("SO:0000448","three_prime_UTR_intron","SO:0000446");
		tree.newTerm("SO:0000449","random_sequence","SO:0000351");
		tree.newTerm("SO:0000450","interband","SO:0000830");
		tree.newTerm("SO:0000451","gene_with_polyadenylated_mRNA","SO:0001217");
		tree.newTerm("SO:0000453","chromosomal_transposition","SO:1000183");
		tree.newTerm("SO:0000454","rasiRNA","SO:0000655");
		tree.newTerm("SO:0000455","gene_with_mRNA_with_frameshift","SO:0001217");
		tree.newTerm("SO:0000456","recombinationally_rearranged_gene","SO:0000704");
		tree.newTerm("SO:0000457","interchromosomal_duplication","SO:1000037");
		tree.newTerm("SO:0000458","D_gene_segment","SO:0000460");
		tree.newTerm("SO:0000459","gene_with_trans_spliced_transcript","SO:0000704");
		tree.newTerm("SO:0000460","vertebrate_immunoglobulin_T_cell_receptor_segment","SO:0000301");
		tree.newTerm("SO:0000461","inversion_derived_bipartite_deficiency","SO:1000029");
		tree.newTerm("SO:0000462","pseudogenic_region","SO:0001411");
		tree.newTerm("SO:0000463","encodes_alternately_spliced_transcripts","SO:0000401");
		tree.newTerm("SO:0000464","decayed_exon","SO:0000462");
		tree.newTerm("SO:0000465","inversion_derived_deficiency_plus_duplication","SO:1000029","SO:1000038");
		tree.newTerm("SO:0000466","V_gene_segment","SO:0000460");
		tree.newTerm("SO:0000467","post_translationally_regulated_by_protein_stability","SO:0000130");
		tree.newTerm("SO:0000468","golden_path_fragment","SO:0000143");
		tree.newTerm("SO:0000469","post_translationally_regulated_by_protein_modification","SO:0000130");
		tree.newTerm("SO:0000470","J_gene_segment","SO:0000460");
		tree.newTerm("SO:0000471","autoregulated","SO:0000123");
		tree.newTerm("SO:0000472","tiling_path","SO:0000353");
		tree.newTerm("SO:0000473","negatively_autoregulated","SO:0000126","SO:0000471");
		tree.newTerm("SO:0000474","tiling_path_fragment","SO:0000143");
		tree.newTerm("SO:0000475","positively_autoregulated","SO:0000125","SO:0000471");
		tree.newTerm("SO:0000476","contig_read","SO:0000150");
		tree.newTerm("SO:0000478","C_gene_segment","SO:0000460");
		tree.newTerm("SO:0000479","trans_spliced_transcript","SO:0000673");
		tree.newTerm("SO:0000480","tiling_path_clone","SO:0000151","SO:0000474");
		tree.newTerm("SO:0000481","terminal_inverted_repeat","SO:0000294");
		tree.newTerm("SO:0000482","vertebrate_immunoglobulin_T_cell_receptor_gene_cluster","SO:0000301");
		tree.newTerm("SO:0000483","nc_primary_transcript","SO:0000185");
		tree.newTerm("SO:0000484","three_prime_coding_exon_noncoding_region","SO:0001214");
		tree.newTerm("SO:0000485","DJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000486","five_prime_coding_exon_noncoding_region","SO:0001214");
		tree.newTerm("SO:0000487","VDJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000488","VDJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000489","VJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000490","VJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000491","VJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000492","D_gene_recombination_feature","SO:0000939");
		tree.newTerm("SO:0000493","three_prime_D_heptamer","SO:0000561");
		tree.newTerm("SO:0000494","three_prime_D_nonamer","SO:0000562");
		tree.newTerm("SO:0000495","three_prime_D_spacer","SO:0000563");
		tree.newTerm("SO:0000496","five_prime_D_heptamer","SO:0000561");
		tree.newTerm("SO:0000497","five_prime_D_nonamer","SO:0000562");
		tree.newTerm("SO:0000498","five_prime_D_spacer","SO:0000563");
		tree.newTerm("SO:0000499","virtual_sequence","SO:0000353");
		tree.newTerm("SO:0000500","Hoogsteen_base_pair","SO:0000028");
		tree.newTerm("SO:0000501","reverse_Hoogsteen_base_pair","SO:0000028");
		tree.newTerm("SO:0000504","D_DJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000505","D_DJ_cluster","SO:0000938");
		tree.newTerm("SO:0000506","D_DJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000507","pseudogenic_exon","SO:0000462");
		tree.newTerm("SO:0000508","D_DJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000509","D_J_C_cluster","SO:0000482");
		tree.newTerm("SO:0000510","VD_gene_segment","SO:0000936");
		tree.newTerm("SO:0000511","J_C_cluster","SO:0000482");
		tree.newTerm("SO:0000512","inversion_derived_deficiency_plus_aneuploid","SO:1000029");
		tree.newTerm("SO:0000513","J_cluster","SO:0000482");
		tree.newTerm("SO:0000514","J_nonamer","SO:0000562");
		tree.newTerm("SO:0000515","J_heptamer","SO:0000561");
		tree.newTerm("SO:0000516","pseudogenic_transcript","SO:0000462");
		tree.newTerm("SO:0000517","J_spacer","SO:0000563");
		tree.newTerm("SO:0000518","V_DJ_cluster","SO:0000938");
		tree.newTerm("SO:0000519","V_DJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000520","V_VDJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000521","V_VDJ_cluster","SO:0000938");
		tree.newTerm("SO:0000522","V_VDJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000523","V_VJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000524","V_VJ_cluster","SO:0000938");
		tree.newTerm("SO:0000525","V_VJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000526","V_cluster","SO:0000482");
		tree.newTerm("SO:0000527","V_D_DJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000528","V_D_DJ_cluster","SO:0000938");
		tree.newTerm("SO:0000529","V_D_DJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000530","V_D_DJ_J_cluster","SO:0000938");
		tree.newTerm("SO:0000531","V_D_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000532","V_D_J_cluster","SO:0000938");
		tree.newTerm("SO:0000533","V_heptamer","SO:0000561");
		tree.newTerm("SO:0000534","V_J_cluster","SO:0000482");
		tree.newTerm("SO:0000535","V_J_C_cluster","SO:0000482");
		tree.newTerm("SO:0000536","V_nonamer","SO:0000562");
		tree.newTerm("SO:0000537","V_spacer","SO:0000563");
		tree.newTerm("SO:0000538","V_gene_recombination_feature","SO:0000939");
		tree.newTerm("SO:0000539","DJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000540","DJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000541","VDJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000542","V_DJ_C_cluster","SO:0000938");
		tree.newTerm("SO:0000544","helitron","SO:0000182");
		tree.newTerm("SO:0000545","recoding_pseudoknot","SO:0000591");
		tree.newTerm("SO:0000546","designed_sequence","SO:0000351");
		tree.newTerm("SO:0000547","inversion_derived_bipartite_duplication","SO:1000038");
		tree.newTerm("SO:0000548","gene_with_edited_transcript","SO:0001217");
		tree.newTerm("SO:0000549","inversion_derived_duplication_plus_aneuploid","SO:1000038");
		tree.newTerm("SO:0000550","aneuploid_chromosome","SO:1000183");
		tree.newTerm("SO:0000551","polyA_signal_sequence","SO:0001679");
		tree.newTerm("SO:0000552","Shine_Dalgarno_sequence","SO:0000139");
		tree.newTerm("SO:0000553","polyA_site","SO:0000699");
		tree.newTerm("SO:0000555","five_prime_clip","SO:0000303");
		tree.newTerm("SO:0000556","five_prime_D_recombination_signal_sequence","SO:0000492");
		tree.newTerm("SO:0000557","three_prime_clip","SO:0000303");
		tree.newTerm("SO:0000558","C_cluster","SO:0000482");
		tree.newTerm("SO:0000559","D_cluster","SO:0000482");
		tree.newTerm("SO:0000560","D_J_cluster","SO:0000482");
		tree.newTerm("SO:0000561","heptamer_of_recombination_feature_of_vertebrate_immune_system_gene","SO:0000939");
		tree.newTerm("SO:0000562","nonamer_of_recombination_feature_of_vertebrate_immune_system_gene","SO:0000939");
		tree.newTerm("SO:0000563","vertebrate_immune_system_gene_recombination_spacer","SO:0000301");
		tree.newTerm("SO:0000564","V_DJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000565","V_VDJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000566","V_VJ_J_C_cluster","SO:0000938");
		tree.newTerm("SO:0000567","inversion_derived_aneuploid_chromosome","SO:0000550");
		tree.newTerm("SO:0000568","bidirectional_promoter","SO:0000167");
		tree.newTerm("SO:0000569","retrotransposed","SO:0000733");
		tree.newTerm("SO:0000570","three_prime_D_recombination_signal_sequence","SO:0000492");
		tree.newTerm("SO:0000571","miRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000572","DJ_gene_segment","SO:0000936");
		tree.newTerm("SO:0000573","rRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000574","VDJ_gene_segment","SO:0000936");
		tree.newTerm("SO:0000575","scRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000576","VJ_gene_segment","SO:0000936");
		tree.newTerm("SO:0000577","centromere","SO:0000628");
		tree.newTerm("SO:0000578","snoRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000579","edited_transcript_feature","SO:0000833");
		tree.newTerm("SO:0000580","methylation_guide_snoRNA_primary_transcript","SO:0000232");
		tree.newTerm("SO:0000581","cap","SO:0001411");
		tree.newTerm("SO:0000582","rRNA_cleavage_snoRNA_primary_transcript","SO:0000232");
		tree.newTerm("SO:0000583","pre_edited_region","SO:0000579");
		tree.newTerm("SO:0000584","tmRNA","SO:0000370");
		tree.newTerm("SO:0000585","C_D_box_snoRNA_encoding","SO:0000578");
		tree.newTerm("SO:0000586","tmRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000587","group_I_intron","SO:0000588");
		tree.newTerm("SO:0000588","autocatalytically_spliced_intron","SO:0000188");
		tree.newTerm("SO:0000589","SRP_RNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000590","SRP_RNA","SO:0000655");
		tree.newTerm("SO:0000591","pseudoknot","SO:0000002");
		tree.newTerm("SO:0000592","H_pseudoknot","SO:0000591");
		tree.newTerm("SO:0000593","C_D_box_snoRNA","SO:0000275");
		tree.newTerm("SO:0000594","H_ACA_box_snoRNA","SO:0000275");
		tree.newTerm("SO:0000595","C_D_box_snoRNA_primary_transcript","SO:0000232");
		tree.newTerm("SO:0000596","H_ACA_box_snoRNA_primary_transcript","SO:0000232");
		tree.newTerm("SO:0000602","guide_RNA","SO:0000655");
		tree.newTerm("SO:0000603","group_II_intron","SO:0000588");
		tree.newTerm("SO:0000604","editing_block","SO:0000579");
		tree.newTerm("SO:0000605","intergenic_region","SO:0001411");
		tree.newTerm("SO:0000606","editing_domain","SO:0000579");
		tree.newTerm("SO:0000607","unedited_region","SO:0000579");
		tree.newTerm("SO:0000608","H_ACA_box_snoRNA_encoding","SO:0000578");
		tree.newTerm("SO:0000609","oligo_U_tail","SO:0001411");
		tree.newTerm("SO:0000610","polyA_sequence","SO:0001411");
		tree.newTerm("SO:0000611","branch_site","SO:0000841");
		tree.newTerm("SO:0000612","polypyrimidine_tract","SO:0000841");
		tree.newTerm("SO:0000613","bacterial_RNApol_promoter","SO:0000752","SO:0001203");
		tree.newTerm("SO:0000614","bacterial_terminator","SO:0000141","SO:0000752");
		tree.newTerm("SO:0000615","terminator_of_type_2_RNApol_III_promoter","SO:0000951");
		tree.newTerm("SO:0000616","transcription_end_site","SO:0000835");
		tree.newTerm("SO:0000617","RNApol_III_promoter_type_1","SO:0000171");
		tree.newTerm("SO:0000618","RNApol_III_promoter_type_2","SO:0000171");
		tree.newTerm("SO:0000619","A_box","SO:0001660");
		tree.newTerm("SO:0000620","B_box","SO:0001660");
		tree.newTerm("SO:0000621","RNApol_III_promoter_type_3","SO:0000171");
		tree.newTerm("SO:0000622","C_box","SO:0001660");
		tree.newTerm("SO:0000623","snRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000624","telomere","SO:0000628");
		tree.newTerm("SO:0000625","silencer","SO:0000727");
		tree.newTerm("SO:0000626","chromosomal_regulatory_element","SO:0000830");
		tree.newTerm("SO:0000627","insulator","SO:0001055");
		tree.newTerm("SO:0000628","chromosomal_structural_element","SO:0000830");
		tree.newTerm("SO:0000629","five_prime_open_reading_frame","SO:0000836");
		tree.newTerm("SO:0000630","upstream_AUG_codon","SO:0000837");
		tree.newTerm("SO:0000631","polycistronic_primary_transcript","SO:0000078","SO:0000185");
		tree.newTerm("SO:0000632","monocistronic_primary_transcript","SO:0000185","SO:0000665");
		tree.newTerm("SO:0000633","monocistronic_mRNA","SO:0000234","SO:0000665");
		tree.newTerm("SO:0000634","polycistronic_mRNA","SO:0000078","SO:0000234");
		tree.newTerm("SO:0000635","mini_exon_donor_RNA","SO:0000185");
		tree.newTerm("SO:0000636","spliced_leader_RNA","SO:0000835");
		tree.newTerm("SO:0000637","engineered_plasmid","SO:0000155","SO:0000804");
		tree.newTerm("SO:0000638","transcribed_spacer_region","SO:0000838");
		tree.newTerm("SO:0000639","internal_transcribed_spacer_region","SO:0000638");
		tree.newTerm("SO:0000640","external_transcribed_spacer_region","SO:0000638");
		tree.newTerm("SO:0000641","tetranucleotide_repeat_microsatellite_feature","SO:0000289");
		tree.newTerm("SO:0000642","SRP_RNA_encoding","SO:0000011");
		tree.newTerm("SO:0000643","minisatellite","SO:0000005");
		tree.newTerm("SO:0000644","antisense_RNA","SO:0000655");
		tree.newTerm("SO:0000645","antisense_primary_transcript","SO:0000185");
		tree.newTerm("SO:0000646","siRNA","SO:0000655");
		tree.newTerm("SO:0000647","miRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0000650","small_subunit_rRNA","SO:0000252");
		tree.newTerm("SO:0000651","large_subunit_rRNA","SO:0000252");
		tree.newTerm("SO:0000652","rRNA_5S","SO:0000651");
		tree.newTerm("SO:0000653","rRNA_28S","SO:0000651");
		tree.newTerm("SO:0000654","maxicircle_gene","SO:0000089");
		tree.newTerm("SO:0000655","ncRNA","SO:0000233");
		tree.newTerm("SO:0000656","stRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000657","repeat_region","SO:0001411");
		tree.newTerm("SO:0000658","dispersed_repeat","SO:0000657");
		tree.newTerm("SO:0000659","tmRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000662","spliceosomal_intron","SO:0000188");
		tree.newTerm("SO:0000663","tRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000664","introgressed_chromosome_region","SO:0000830");
		tree.newTerm("SO:0000665","monocistronic_transcript","SO:0000673");
		tree.newTerm("SO:0000666","mobile_intron","SO:0000188","SO:0001037");
		tree.newTerm("SO:0000667","insertion","SO:0001059","SO:0001411");
		tree.newTerm("SO:0000668","EST_match","SO:0000102");
		tree.newTerm("SO:0000669","sequence_rearrangement_feature","SO:0000298");
		tree.newTerm("SO:0000670","chromosome_breakage_sequence","SO:0000669");
		tree.newTerm("SO:0000671","internal_eliminated_sequence","SO:0000669");
		tree.newTerm("SO:0000672","macronucleus_destined_segment","SO:0000669");
		tree.newTerm("SO:0000673","transcript","SO:0000831");
		tree.newTerm("SO:0000676","canonical_three_prime_splice_site","SO:0000164");
		tree.newTerm("SO:0000677","canonical_five_prime_splice_site","SO:0000163");
		tree.newTerm("SO:0000678","non_canonical_three_prime_splice_site","SO:0000164");
		tree.newTerm("SO:0000679","non_canonical_five_prime_splice_site","SO:0000163");
		tree.newTerm("SO:0000680","non_canonical_start_codon","SO:0000318");
		tree.newTerm("SO:0000681","aberrant_processed_transcript","SO:0000673");
		tree.newTerm("SO:0000683","exonic_splice_enhancer","SO:0000344");
		tree.newTerm("SO:0000684","nuclease_sensitive_site","SO:0000059");
		tree.newTerm("SO:0000685","DNAseI_hypersensitive_site","SO:0000322");
		tree.newTerm("SO:0000686","translocation_element","SO:1000044");
		tree.newTerm("SO:0000687","deletion_junction","SO:0000699");
		tree.newTerm("SO:0000688","golden_path","SO:0000353");
		tree.newTerm("SO:0000689","cDNA_match","SO:0000102");
		tree.newTerm("SO:0000690","gene_with_polycistronic_transcript","SO:0000704");
		tree.newTerm("SO:0000691","cleaved_initiator_methionine","SO:0100011");
		tree.newTerm("SO:0000692","gene_with_dicistronic_transcript","SO:0000690");
		tree.newTerm("SO:0000693","gene_with_recoded_mRNA","SO:0001217");
		tree.newTerm("SO:0000694","SNP","SO:0001483");
		tree.newTerm("SO:0000695","reagent","SO:0001409");
		tree.newTerm("SO:0000696","oligo","SO:0000695");
		tree.newTerm("SO:0000697","gene_with_stop_codon_read_through","SO:0000693");
		tree.newTerm("SO:0000698","gene_with_stop_codon_redefined_as_pyrrolysine","SO:0000697");
		tree.newTerm("SO:0000699","junction","SO:0000110");
		tree.newTerm("SO:0000700","remark","SO:0001410");
		tree.newTerm("SO:0000701","possible_base_call_error","SO:0000413");
		tree.newTerm("SO:0000702","possible_assembly_error","SO:0000413");
		tree.newTerm("SO:0000703","experimental_result_region","SO:0000700");
		tree.newTerm("SO:0000704","gene","SO:0001411");
		tree.newTerm("SO:0000705","tandem_repeat","SO:0000657");
		tree.newTerm("SO:0000706","trans_splice_acceptor_site","SO:0001420");
		tree.newTerm("SO:0000707","trans_splice_donor_site","SO:0001420");
		tree.newTerm("SO:0000708","SL1_acceptor_site","SO:0000706");
		tree.newTerm("SO:0000709","SL2_acceptor_site","SO:0000706");
		tree.newTerm("SO:0000710","gene_with_stop_codon_redefined_as_selenocysteine","SO:0000697");
		tree.newTerm("SO:0000711","gene_with_mRNA_recoded_by_translational_bypass","SO:0000693");
		tree.newTerm("SO:0000712","gene_with_transcript_with_translational_frameshift","SO:0000693");
		tree.newTerm("SO:0000713","DNA_motif","SO:0000714");
		tree.newTerm("SO:0000714","nucleotide_motif","SO:0001683");
		tree.newTerm("SO:0000715","RNA_motif","SO:0000714");
		tree.newTerm("SO:0000716","dicistronic_mRNA","SO:0000079","SO:0000634");
		tree.newTerm("SO:0000717","reading_frame","SO:0001410");
		tree.newTerm("SO:0000718","blocked_reading_frame","SO:0000717");
		tree.newTerm("SO:0000719","ultracontig","SO:0001876");
		tree.newTerm("SO:0000720","foreign_transposable_element","SO:0000101");
		tree.newTerm("SO:0000721","gene_with_dicistronic_primary_transcript","SO:0000692");
		tree.newTerm("SO:0000722","gene_with_dicistronic_mRNA","SO:0000692");
		tree.newTerm("SO:0000723","iDNA","SO:0000298");
		tree.newTerm("SO:0000724","oriT","SO:0000296");
		tree.newTerm("SO:0000725","transit_peptide","SO:0001527");
		tree.newTerm("SO:0000726","repeat_unit","SO:0001411");
		tree.newTerm("SO:0000727","CRM","SO:0001055");
		tree.newTerm("SO:0000728","intein","SO:0100011");
		tree.newTerm("SO:0000729","intein_containing","SO:0000010");
		tree.newTerm("SO:0000730","gap","SO:0000143");
		tree.newTerm("SO:0000731","fragmentary","SO:0000905");
		tree.newTerm("SO:0000732","predicted","SO:0000905");
		tree.newTerm("SO:0000733","feature_attribute","SO:0000400");
		tree.newTerm("SO:0000734","exemplar_mRNA","SO:0000234");
		tree.newTerm("SO:0000735","sequence_location","SO:0000400");
		tree.newTerm("SO:0000736","organelle_sequence","SO:0000735");
		tree.newTerm("SO:0000737","mitochondrial_sequence","SO:0000736");
		tree.newTerm("SO:0000738","nuclear_sequence","SO:0000736");
		tree.newTerm("SO:0000739","nucleomorphic_sequence","SO:0000736");
		tree.newTerm("SO:0000740","plastid_sequence","SO:0000736");
		tree.newTerm("SO:0000741","kinetoplast","SO:0001026");
		tree.newTerm("SO:0000742","maxicircle","SO:0001235");
		tree.newTerm("SO:0000743","apicoplast_sequence","SO:0000740");
		tree.newTerm("SO:0000744","chromoplast_sequence","SO:0000740");
		tree.newTerm("SO:0000745","chloroplast_sequence","SO:0000740");
		tree.newTerm("SO:0000746","cyanelle_sequence","SO:0000740");
		tree.newTerm("SO:0000747","leucoplast_sequence","SO:0000740");
		tree.newTerm("SO:0000748","proplastid_sequence","SO:0000740");
		tree.newTerm("SO:0000749","plasmid_location","SO:0000735");
		tree.newTerm("SO:0000750","amplification_origin","SO:0000296");
		tree.newTerm("SO:0000751","proviral_location","SO:0000735");
		tree.newTerm("SO:0000752","gene_group_regulatory_region","SO:0001679");
		tree.newTerm("SO:0000753","clone_insert","SO:0000695");
		tree.newTerm("SO:0000754","lambda_vector","SO:0000440");
		tree.newTerm("SO:0000755","plasmid_vector","SO:0000440");
		tree.newTerm("SO:0000756","cDNA","SO:0000352");
		tree.newTerm("SO:0000757","single_stranded_cDNA","SO:0000756");
		tree.newTerm("SO:0000758","double_stranded_cDNA","SO:0000756");
		tree.newTerm("SO:0000766","pyrrolysyl_tRNA","SO:0000253");
		tree.newTerm("SO:0000768","episome","SO:0000155");
		tree.newTerm("SO:0000769","tmRNA_coding_piece","SO:0000847");
		tree.newTerm("SO:0000770","tmRNA_acceptor_piece","SO:0000847");
		tree.newTerm("SO:0000771","QTL","SO:0001411");
		tree.newTerm("SO:0000772","genomic_island","SO:0001039");
		tree.newTerm("SO:0000773","pathogenic_island","SO:0000772");
		tree.newTerm("SO:0000774","metabolic_island","SO:0000772");
		tree.newTerm("SO:0000775","adaptive_island","SO:0000772");
		tree.newTerm("SO:0000776","symbiosis_island","SO:0000772");
		tree.newTerm("SO:0000777","pseudogenic_rRNA","SO:0000516");
		tree.newTerm("SO:0000778","pseudogenic_tRNA","SO:0000516");
		tree.newTerm("SO:0000779","engineered_episome","SO:0000637","SO:0000768");
		tree.newTerm("SO:0000781","transgenic","SO:0000733");
		tree.newTerm("SO:0000782","natural","SO:0000733");
		tree.newTerm("SO:0000783","engineered","SO:0000733");
		tree.newTerm("SO:0000784","foreign","SO:0000733");
		tree.newTerm("SO:0000785","cloned_region","SO:0000695");
		tree.newTerm("SO:0000789","validated","SO:0000905");
		tree.newTerm("SO:0000790","invalidated","SO:0000905");
		tree.newTerm("SO:0000794","engineered_rescue_region","SO:0000411","SO:0000804");
		tree.newTerm("SO:0000795","rescue_mini_gene","SO:0000815");
		tree.newTerm("SO:0000796","transgenic_transposable_element","SO:0000101");
		tree.newTerm("SO:0000797","natural_transposable_element","SO:0000101","SO:0001038");
		tree.newTerm("SO:0000798","engineered_transposable_element","SO:0000101","SO:0000804");
		tree.newTerm("SO:0000799","engineered_foreign_transposable_element","SO:0000720","SO:0000798","SO:0000805");
		tree.newTerm("SO:0000800","assortment_derived_duplication","SO:0001504");
		tree.newTerm("SO:0000801","assortment_derived_deficiency_plus_duplication","SO:0001504");
		tree.newTerm("SO:0000802","assortment_derived_deficiency","SO:0001504");
		tree.newTerm("SO:0000803","assortment_derived_aneuploid","SO:0001504");
		tree.newTerm("SO:0000804","engineered_region","SO:0001409");
		tree.newTerm("SO:0000805","engineered_foreign_region","SO:0000804");
		tree.newTerm("SO:0000806","fusion","SO:0000733");
		tree.newTerm("SO:0000807","engineered_tag","SO:0000324","SO:0000804");
		tree.newTerm("SO:0000808","validated_cDNA_clone","SO:0000317");
		tree.newTerm("SO:0000809","invalidated_cDNA_clone","SO:0000317");
		tree.newTerm("SO:0000810","chimeric_cDNA_clone","SO:0000809");
		tree.newTerm("SO:0000811","genomically_contaminated_cDNA_clone","SO:0000809");
		tree.newTerm("SO:0000812","polyA_primed_cDNA_clone","SO:0000809");
		tree.newTerm("SO:0000813","partially_processed_cDNA_clone","SO:0000809");
		tree.newTerm("SO:0000814","rescue","SO:0000733");
		tree.newTerm("SO:0000815","mini_gene","SO:0000236");
		tree.newTerm("SO:0000816","rescue_gene","SO:0000704");
		tree.newTerm("SO:0000817","wild_type","SO:0000733");
		tree.newTerm("SO:0000818","wild_type_rescue_gene","SO:0000816");
		tree.newTerm("SO:0000819","mitochondrial_chromosome","SO:0000340");
		tree.newTerm("SO:0000820","chloroplast_chromosome","SO:0000340");
		tree.newTerm("SO:0000821","chromoplast_chromosome","SO:0000340");
		tree.newTerm("SO:0000822","cyanelle_chromosome","SO:0000340");
		tree.newTerm("SO:0000823","leucoplast_chromosome","SO:0000340");
		tree.newTerm("SO:0000824","macronuclear_chromosome","SO:0000340");
		tree.newTerm("SO:0000825","micronuclear_chromosome","SO:0000340");
		tree.newTerm("SO:0000828","nuclear_chromosome","SO:0000340");
		tree.newTerm("SO:0000829","nucleomorphic_chromosome","SO:0000340");
		tree.newTerm("SO:0000830","chromosome_part","SO:0001411");
		tree.newTerm("SO:0000831","gene_member_region","SO:0001411");
		tree.newTerm("SO:0000833","transcript_region","SO:0001411");
		tree.newTerm("SO:0000834","mature_transcript_region","SO:0000833");
		tree.newTerm("SO:0000835","primary_transcript_region","SO:0000833");
		tree.newTerm("SO:0000836","mRNA_region","SO:0000834");
		tree.newTerm("SO:0000837","UTR_region","SO:0000836");
		tree.newTerm("SO:0000838","rRNA_primary_transcript_region","SO:0000835");
		tree.newTerm("SO:0000839","polypeptide_region","SO:0001411");
		tree.newTerm("SO:0000840","repeat_component","SO:0001412");
		tree.newTerm("SO:0000841","spliceosomal_intron_region","SO:0000835");
		tree.newTerm("SO:0000842","gene_component_region","SO:0001411");
		tree.newTerm("SO:0000847","tmRNA_region","SO:0000834");
		tree.newTerm("SO:0000848","LTR_component","SO:0000840");
		tree.newTerm("SO:0000849","three_prime_LTR_component","SO:0000848");
		tree.newTerm("SO:0000850","five_prime_LTR_component","SO:0000848");
		tree.newTerm("SO:0000851","CDS_region","SO:0000836");
		tree.newTerm("SO:0000852","exon_region","SO:0000833");
		tree.newTerm("SO:0000853","homologous_region","SO:0000330");
		tree.newTerm("SO:0000854","paralogous_region","SO:0000853");
		tree.newTerm("SO:0000855","orthologous_region","SO:0000853");
		tree.newTerm("SO:0000856","conserved","SO:0000733");
		tree.newTerm("SO:0000857","homologous","SO:0000856");
		tree.newTerm("SO:0000858","orthologous","SO:0000857");
		tree.newTerm("SO:0000859","paralogous","SO:0000857");
		tree.newTerm("SO:0000860","syntenic","SO:0000856");
		tree.newTerm("SO:0000861","capped_primary_transcript","SO:0000185");
		tree.newTerm("SO:0000862","capped_mRNA","SO:0000234");
		tree.newTerm("SO:0000863","mRNA_attribute","SO:0000237");
		tree.newTerm("SO:0000864","exemplar","SO:0000863");
		tree.newTerm("SO:0000865","frameshift","SO:0000863");
		tree.newTerm("SO:0000866","minus_1_frameshift","SO:0000865");
		tree.newTerm("SO:0000867","minus_2_frameshift","SO:0000865");
		tree.newTerm("SO:0000868","plus_1_frameshift","SO:0000865");
		tree.newTerm("SO:0000869","plus_2_framshift","SO:0000865");
		tree.newTerm("SO:0000870","trans_spliced","SO:0000237");
		tree.newTerm("SO:0000871","polyadenylated_mRNA","SO:0000234");
		tree.newTerm("SO:0000872","trans_spliced_mRNA","SO:0000234","SO:0000479");
		tree.newTerm("SO:0000873","edited_transcript","SO:0000673");
		tree.newTerm("SO:0000874","edited_transcript_by_A_to_I_substitution","SO:0000873");
		tree.newTerm("SO:0000875","bound_by_protein","SO:0000277");
		tree.newTerm("SO:0000876","bound_by_nucleic_acid","SO:0000277");
		tree.newTerm("SO:0000877","alternatively_spliced","SO:0000237");
		tree.newTerm("SO:0000878","monocistronic","SO:0000237");
		tree.newTerm("SO:0000879","dicistronic","SO:0000880");
		tree.newTerm("SO:0000880","polycistronic","SO:0000237");
		tree.newTerm("SO:0000881","recoded","SO:0000863");
		tree.newTerm("SO:0000882","codon_redefined","SO:0000881");
		tree.newTerm("SO:0000883","stop_codon_read_through","SO:0000145");
		tree.newTerm("SO:0000884","stop_codon_redefined_as_pyrrolysine","SO:0000883");
		tree.newTerm("SO:0000885","stop_codon_redefined_as_selenocysteine","SO:0000883");
		tree.newTerm("SO:0000886","recoded_by_translational_bypass","SO:0000881");
		tree.newTerm("SO:0000887","translationally_frameshifted","SO:0000881");
		tree.newTerm("SO:0000888","maternally_imprinted_gene","SO:0000898");
		tree.newTerm("SO:0000889","paternally_imprinted_gene","SO:0000898");
		tree.newTerm("SO:0000890","post_translationally_regulated_gene","SO:0000704");
		tree.newTerm("SO:0000891","negatively_autoregulated_gene","SO:0000704");
		tree.newTerm("SO:0000892","positively_autoregulated_gene","SO:0000704");
		tree.newTerm("SO:0000893","silenced","SO:0000126");
		tree.newTerm("SO:0000894","silenced_by_DNA_modification","SO:0000893");
		tree.newTerm("SO:0000895","silenced_by_DNA_methylation","SO:0000894");
		tree.newTerm("SO:0000896","translationally_regulated_gene","SO:0000704");
		tree.newTerm("SO:0000897","allelically_excluded_gene","SO:0000898");
		tree.newTerm("SO:0000898","epigenetically_modified_gene","SO:0000704","SO:0001720");
		tree.newTerm("SO:0000902","transgene","SO:0000704");
		tree.newTerm("SO:0000903","endogenous_retroviral_sequence","SO:0000751");
		tree.newTerm("SO:0000904","rearranged_at_DNA_level","SO:0000133");
		tree.newTerm("SO:0000905","status","SO:0000733");
		tree.newTerm("SO:0000906","independently_known","SO:0000905");
		tree.newTerm("SO:0000907","supported_by_sequence_similarity","SO:0000732");
		tree.newTerm("SO:0000908","supported_by_domain_match","SO:0000907");
		tree.newTerm("SO:0000909","supported_by_EST_or_cDNA","SO:0000907");
		tree.newTerm("SO:0000910","orphan","SO:0000732");
		tree.newTerm("SO:0000911","predicted_by_ab_initio_computation","SO:0000732");
		tree.newTerm("SO:0000912","asx_turn","SO:0001128");
		tree.newTerm("SO:0000913","cloned_cDNA_insert","SO:0000753");
		tree.newTerm("SO:0000914","cloned_genomic_insert","SO:0000753");
		tree.newTerm("SO:0000915","engineered_insert","SO:0000753","SO:0000804");
		tree.newTerm("SO:0000929","edited_mRNA","SO:0000234","SO:0000873");
		tree.newTerm("SO:0000930","guide_RNA_region","SO:0000834");
		tree.newTerm("SO:0000931","anchor_region","SO:0000930");
		tree.newTerm("SO:0000932","pre_edited_mRNA","SO:0000120");
		tree.newTerm("SO:0000933","intermediate","SO:0000733");
		tree.newTerm("SO:0000934","miRNA_target_site","SO:0001655");
		tree.newTerm("SO:0000935","edited_CDS","SO:0000316");
		tree.newTerm("SO:0000936","vertebrate_immunoglobulin_T_cell_receptor_rearranged_segment","SO:0000301");
		tree.newTerm("SO:0000938","vertebrate_immunoglobulin_T_cell_receptor_rearranged_gene_cluster","SO:0000301");
		tree.newTerm("SO:0000939","vertebrate_immune_system_gene_recombination_signal_feature","SO:0000301");
		tree.newTerm("SO:0000940","recombinationally_rearranged","SO:0000733");
		tree.newTerm("SO:0000941","recombinationally_rearranged_vertebrate_immune_system_gene","SO:0000456");
		tree.newTerm("SO:0000942","attP_site","SO:0000946");
		tree.newTerm("SO:0000943","attB_site","SO:0000946");
		tree.newTerm("SO:0000944","attL_site","SO:0000946");
		tree.newTerm("SO:0000945","attR_site","SO:0000946");
		tree.newTerm("SO:0000946","integration_excision_site","SO:0000342");
		tree.newTerm("SO:0000947","resolution_site","SO:0000342");
		tree.newTerm("SO:0000948","inversion_site","SO:0000342");
		tree.newTerm("SO:0000949","dif_site","SO:0000947");
		tree.newTerm("SO:0000950","attC_site","SO:0000946");
		tree.newTerm("SO:0000951","eukaryotic_terminator","SO:0000141");
		tree.newTerm("SO:0000952","oriV","SO:0000296");
		tree.newTerm("SO:0000953","oriC","SO:0000296");
		tree.newTerm("SO:0000954","DNA_chromosome","SO:0000340");
		tree.newTerm("SO:0000955","double_stranded_DNA_chromosome","SO:0000954");
		tree.newTerm("SO:0000956","single_stranded_DNA_chromosome","SO:0000954");
		tree.newTerm("SO:0000957","linear_double_stranded_DNA_chromosome","SO:0000955");
		tree.newTerm("SO:0000958","circular_double_stranded_DNA_chromosome","SO:0000955");
		tree.newTerm("SO:0000959","linear_single_stranded_DNA_chromosome","SO:0000956");
		tree.newTerm("SO:0000960","circular_single_stranded_DNA_chromosome","SO:0000956");
		tree.newTerm("SO:0000961","RNA_chromosome","SO:0000340");
		tree.newTerm("SO:0000962","single_stranded_RNA_chromosome","SO:0000961");
		tree.newTerm("SO:0000963","linear_single_stranded_RNA_chromosome","SO:0000962");
		tree.newTerm("SO:0000964","linear_double_stranded_RNA_chromosome","SO:0000965");
		tree.newTerm("SO:0000965","double_stranded_RNA_chromosome","SO:0000961");
		tree.newTerm("SO:0000966","circular_single_stranded_RNA_chromosome","SO:0000962");
		tree.newTerm("SO:0000967","circular_double_stranded_RNA_chromosome","SO:0000965");
		tree.newTerm("SO:0000973","insertion_sequence","SO:0000208");
		tree.newTerm("SO:0000975","minicircle_gene","SO:0000089");
		tree.newTerm("SO:0000976","cryptic","SO:0000733");
		tree.newTerm("SO:0000977","anchor_binding_site","SO:0000833");
		tree.newTerm("SO:0000978","template_region","SO:0000930");
		tree.newTerm("SO:0000979","gRNA_encoding","SO:0000011");
		tree.newTerm("SO:0000980","minicircle","SO:0001235");
		tree.newTerm("SO:0000981","rho_dependent_bacterial_terminator","SO:0000614");
		tree.newTerm("SO:0000982","rho_independent_bacterial_terminator","SO:0000614");
		tree.newTerm("SO:0000983","strand_attribute","SO:0000733");
		tree.newTerm("SO:0000984","single","SO:0000983");
		tree.newTerm("SO:0000985","double","SO:0000983");
		tree.newTerm("SO:0000986","topology_attribute","SO:0000443");
		tree.newTerm("SO:0000987","linear","SO:0000986");
		tree.newTerm("SO:0000988","circular","SO:0000986");
		tree.newTerm("SO:0000989","class_II_RNA","SO:0000655");
		tree.newTerm("SO:0000990","class_I_RNA","SO:0000655");
		tree.newTerm("SO:0000991","genomic_DNA","SO:0000352");
		tree.newTerm("SO:0000992","BAC_cloned_genomic_insert","SO:0000914");
		tree.newTerm("SO:0000993","consensus","SO:0000905");
		tree.newTerm("SO:0000994","consensus_region","SO:0001410");
		tree.newTerm("SO:0000995","consensus_mRNA","SO:0000234","SO:0000994");
		tree.newTerm("SO:0000996","predicted_gene","SO:0000704");
		tree.newTerm("SO:0000997","gene_fragment","SO:0000842");
		tree.newTerm("SO:0000998","recursive_splice_site","SO:0001419");
		tree.newTerm("SO:0000999","BAC_end","SO:0000150");
		tree.newTerm("SO:0001000","rRNA_16S","SO:0000650");
		tree.newTerm("SO:0001001","rRNA_23S","SO:0000651");
		tree.newTerm("SO:0001002","rRNA_25S","SO:0000651");
		tree.newTerm("SO:0001003","solo_LTR","SO:0000286");
		tree.newTerm("SO:0001004","low_complexity","SO:0000905");
		tree.newTerm("SO:0001005","low_complexity_region","SO:0001410");
		tree.newTerm("SO:0001006","prophage","SO:0000113");
		tree.newTerm("SO:0001007","cryptic_prophage","SO:0000772");
		tree.newTerm("SO:0001008","tetraloop","SO:0000313");
		tree.newTerm("SO:0001009","DNA_constraint_sequence","SO:0000442");
		tree.newTerm("SO:0001010","i_motif","SO:0000142");
		tree.newTerm("SO:0001011","PNA_oligo","SO:0001247");
		tree.newTerm("SO:0001012","DNAzyme","SO:0000696");
		tree.newTerm("SO:0001013","MNP","SO:0002007");
		tree.newTerm("SO:0001014","intron_domain","SO:0000835");
		tree.newTerm("SO:0001015","wobble_base_pair","SO:0000028");
		tree.newTerm("SO:0001016","internal_guide_sequence","SO:0001014");
		tree.newTerm("SO:0001017","silent_mutation","SO:0001878");
		tree.newTerm("SO:0001018","epitope","SO:0000409");
		tree.newTerm("SO:0001019","copy_number_variation","SO:0000248");
		tree.newTerm("SO:0001021","chromosome_breakpoint","SO:0000699");
		tree.newTerm("SO:0001022","inversion_breakpoint","SO:0001021");
		tree.newTerm("SO:0001023","allele","SO:0001507");
		tree.newTerm("SO:0001024","haplotype","SO:0001507");
		tree.newTerm("SO:0001025","polymorphic_sequence_variant","SO:0001023");
		tree.newTerm("SO:0001026","genome","SO:0001260");
		tree.newTerm("SO:0001027","genotype","SO:0001507");
		tree.newTerm("SO:0001028","diplotype","SO:0001507");
		tree.newTerm("SO:0001029","direction_attribute","SO:0000733");
		tree.newTerm("SO:0001030","forward","SO:0001029");
		tree.newTerm("SO:0001031","reverse","SO:0001029");
		tree.newTerm("SO:0001032","mitochondrial_DNA","SO:0000737");
		tree.newTerm("SO:0001033","chloroplast_DNA","SO:0000745");
		tree.newTerm("SO:0001034","miRtron","SO:0001014");
		tree.newTerm("SO:0001035","piRNA","SO:0000655");
		tree.newTerm("SO:0001036","arginyl_tRNA","SO:0000253");
		tree.newTerm("SO:0001037","mobile_genetic_element","SO:0001411");
		tree.newTerm("SO:0001038","extrachromosomal_mobile_genetic_element","SO:0001037");
		tree.newTerm("SO:0001039","integrated_mobile_genetic_element","SO:0001037");
		tree.newTerm("SO:0001040","integrated_plasmid","SO:0001039");
		tree.newTerm("SO:0001041","viral_sequence","SO:0001038","SO:0001235");
		tree.newTerm("SO:0001042","phage_sequence","SO:0001041");
		tree.newTerm("SO:0001043","attCtn_site","SO:0000946");
		tree.newTerm("SO:0001044","nuclear_mt_pseudogene","SO:0001760");
		tree.newTerm("SO:0001045","cointegrated_plasmid","SO:0001039");
		tree.newTerm("SO:0001046","IRLinv_site","SO:0001048");
		tree.newTerm("SO:0001047","IRRinv_site","SO:0001048");
		tree.newTerm("SO:0001048","inversion_site_part","SO:0000342");
		tree.newTerm("SO:0001049","defective_conjugative_transposon","SO:0000772");
		tree.newTerm("SO:0001050","repeat_fragment","SO:0000840");
		tree.newTerm("SO:0001054","transposon_fragment","SO:0000840");
		tree.newTerm("SO:0001055","transcriptional_cis_regulatory_region","SO:0001679");
		tree.newTerm("SO:0001056","splicing_regulatory_region","SO:0001679");
		tree.newTerm("SO:0001058","promoter_targeting_sequence","SO:0001055");
		tree.newTerm("SO:0001059","sequence_alteration","SO:0002072");
		tree.newTerm("SO:0001060","sequence_variant");
		tree.newTerm("SO:0001061","propeptide_cleavage_site","SO:0100011");
		tree.newTerm("SO:0001062","propeptide","SO:0100011");
		tree.newTerm("SO:0001063","immature_peptide_region","SO:0000839");
		tree.newTerm("SO:0001064","active_peptide","SO:0000419");
		tree.newTerm("SO:0001066","compositionally_biased_region_of_peptide","SO:0000839");
		tree.newTerm("SO:0001067","polypeptide_motif","SO:0100021");
		tree.newTerm("SO:0001068","polypeptide_repeat","SO:0100021");
		tree.newTerm("SO:0001070","polypeptide_structural_region","SO:0000839");
		tree.newTerm("SO:0001071","membrane_structure","SO:0001070");
		tree.newTerm("SO:0001072","extramembrane_polypeptide_region","SO:0001070");
		tree.newTerm("SO:0001073","cytoplasmic_polypeptide_region","SO:0001072");
		tree.newTerm("SO:0001074","non_cytoplasmic_polypeptide_region","SO:0001072");
		tree.newTerm("SO:0001075","intramembrane_polypeptide_region","SO:0001070");
		tree.newTerm("SO:0001076","membrane_peptide_loop","SO:0001075");
		tree.newTerm("SO:0001077","transmembrane_polypeptide_region","SO:0001075");
		tree.newTerm("SO:0001078","polypeptide_secondary_structure","SO:0001070");
		tree.newTerm("SO:0001079","polypeptide_structural_motif","SO:0001070");
		tree.newTerm("SO:0001080","coiled_coil","SO:0001079");
		tree.newTerm("SO:0001081","helix_turn_helix","SO:0001079");
		tree.newTerm("SO:0001082","polypeptide_sequencing_information","SO:0000700");
		tree.newTerm("SO:0001083","non_adjacent_residues","SO:0001082");
		tree.newTerm("SO:0001084","non_terminal_residue","SO:0001082");
		tree.newTerm("SO:0001085","sequence_conflict","SO:0001082");
		tree.newTerm("SO:0001086","sequence_uncertainty","SO:0001082");
		tree.newTerm("SO:0001089","post_translationally_modified_region","SO:0100001");
		tree.newTerm("SO:0001092","polypeptide_metal_contact","SO:0001656","SO:0100002");
		tree.newTerm("SO:0001093","protein_protein_contact","SO:0000410","SO:0100002");
		tree.newTerm("SO:0001094","polypeptide_calcium_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001095","polypeptide_cobalt_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001096","polypeptide_copper_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001097","polypeptide_iron_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001098","polypeptide_magnesium_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001099","polypeptide_manganese_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001100","polypeptide_molybdenum_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001101","polypeptide_nickel_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001102","polypeptide_tungsten_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001103","polypeptide_zinc_ion_contact_site","SO:0001092");
		tree.newTerm("SO:0001104","catalytic_residue","SO:0001237");
		tree.newTerm("SO:0001105","polypeptide_ligand_contact","SO:0001657","SO:0100002");
		tree.newTerm("SO:0001106","asx_motif","SO:0001078");
		tree.newTerm("SO:0001107","beta_bulge","SO:0001078");
		tree.newTerm("SO:0001108","beta_bulge_loop","SO:0001078");
		tree.newTerm("SO:0001109","beta_bulge_loop_five","SO:0001108");
		tree.newTerm("SO:0001110","beta_bulge_loop_six","SO:0001108");
		tree.newTerm("SO:0001111","beta_strand","SO:0001078");
		tree.newTerm("SO:0001112","antiparallel_beta_strand","SO:0001111");
		tree.newTerm("SO:0001113","parallel_beta_strand","SO:0001111");
		tree.newTerm("SO:0001114","peptide_helix","SO:0001078");
		tree.newTerm("SO:0001115","left_handed_peptide_helix","SO:0001114");
		tree.newTerm("SO:0001116","right_handed_peptide_helix","SO:0001114");
		tree.newTerm("SO:0001117","alpha_helix","SO:0001116");
		tree.newTerm("SO:0001118","pi_helix","SO:0001116");
		tree.newTerm("SO:0001119","three_ten_helix","SO:0001116");
		tree.newTerm("SO:0001120","polypeptide_nest_motif","SO:0001078");
		tree.newTerm("SO:0001121","polypeptide_nest_left_right_motif","SO:0001120");
		tree.newTerm("SO:0001122","polypeptide_nest_right_left_motif","SO:0001120");
		tree.newTerm("SO:0001123","schellmann_loop","SO:0001078");
		tree.newTerm("SO:0001124","schellmann_loop_seven","SO:0001123");
		tree.newTerm("SO:0001125","schellmann_loop_six","SO:0001123");
		tree.newTerm("SO:0001126","serine_threonine_motif","SO:0001078");
		tree.newTerm("SO:0001127","serine_threonine_staple_motif","SO:0001078");
		tree.newTerm("SO:0001128","polypeptide_turn_motif","SO:0001078");
		tree.newTerm("SO:0001129","asx_turn_left_handed_type_one","SO:0000912");
		tree.newTerm("SO:0001130","asx_turn_left_handed_type_two","SO:0000912");
		tree.newTerm("SO:0001131","asx_turn_right_handed_type_two","SO:0000912");
		tree.newTerm("SO:0001132","asx_turn_right_handed_type_one","SO:0000912");
		tree.newTerm("SO:0001133","beta_turn","SO:0001128");
		tree.newTerm("SO:0001134","beta_turn_left_handed_type_one","SO:0001133");
		tree.newTerm("SO:0001135","beta_turn_left_handed_type_two","SO:0001133");
		tree.newTerm("SO:0001136","beta_turn_right_handed_type_one","SO:0001133");
		tree.newTerm("SO:0001137","beta_turn_right_handed_type_two","SO:0001133");
		tree.newTerm("SO:0001138","gamma_turn","SO:0001128");
		tree.newTerm("SO:0001139","gamma_turn_classic","SO:0001138");
		tree.newTerm("SO:0001140","gamma_turn_inverse","SO:0001138");
		tree.newTerm("SO:0001141","serine_threonine_turn","SO:0001128");
		tree.newTerm("SO:0001142","st_turn_left_handed_type_one","SO:0001141");
		tree.newTerm("SO:0001143","st_turn_left_handed_type_two","SO:0001141");
		tree.newTerm("SO:0001144","st_turn_right_handed_type_one","SO:0001141");
		tree.newTerm("SO:0001145","st_turn_right_handed_type_two","SO:0001141");
		tree.newTerm("SO:0001146","polypeptide_variation_site","SO:0000839");
		tree.newTerm("SO:0001147","natural_variant_site","SO:0001146");
		tree.newTerm("SO:0001148","mutated_variant_site","SO:0001146");
		tree.newTerm("SO:0001149","alternate_sequence_site","SO:0001146");
		tree.newTerm("SO:0001150","beta_turn_type_six","SO:0001133");
		tree.newTerm("SO:0001151","beta_turn_type_six_a","SO:0001150");
		tree.newTerm("SO:0001152","beta_turn_type_six_a_one","SO:0001151");
		tree.newTerm("SO:0001153","beta_turn_type_six_a_two","SO:0001151");
		tree.newTerm("SO:0001154","beta_turn_type_six_b","SO:0001150");
		tree.newTerm("SO:0001155","beta_turn_type_eight","SO:0001133");
		tree.newTerm("SO:0001156","DRE_motif","SO:0000713");
		tree.newTerm("SO:0001157","DMv4_motif","SO:0001659");
		tree.newTerm("SO:0001158","E_box_motif","SO:0000713");
		tree.newTerm("SO:0001159","DMv5_motif","SO:0001659");
		tree.newTerm("SO:0001160","DMv3_motif","SO:0001659");
		tree.newTerm("SO:0001161","DMv2_motif","SO:0001659");
		tree.newTerm("SO:0001162","MTE","SO:0001660");
		tree.newTerm("SO:0001163","INR1_motif","SO:0000713");
		tree.newTerm("SO:0001164","DPE1_motif","SO:0001659");
		tree.newTerm("SO:0001165","DMv1_motif","SO:0001659");
		tree.newTerm("SO:0001166","GAGA_motif","SO:0000713");
		tree.newTerm("SO:0001167","NDM2_motif","SO:0001659");
		tree.newTerm("SO:0001168","NDM3_motif","SO:0001659");
		tree.newTerm("SO:0001169","ds_RNA_viral_sequence","SO:0001041");
		tree.newTerm("SO:0001170","polinton","SO:0000208");
		tree.newTerm("SO:0001171","rRNA_21S","SO:0000651");
		tree.newTerm("SO:0001172","tRNA_region","SO:0000834");
		tree.newTerm("SO:0001173","anticodon_loop","SO:0001172");
		tree.newTerm("SO:0001174","anticodon","SO:0001172");
		tree.newTerm("SO:0001175","CCA_tail","SO:0001172");
		tree.newTerm("SO:0001176","DHU_loop","SO:0001172");
		tree.newTerm("SO:0001177","T_loop","SO:0001172");
		tree.newTerm("SO:0001178","pyrrolysine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0001179","U3_snoRNA","SO:0000593");
		tree.newTerm("SO:0001180","AU_rich_element","SO:0000837");
		tree.newTerm("SO:0001181","Bruno_response_element","SO:0000837");
		tree.newTerm("SO:0001182","iron_responsive_element","SO:0000837");
		tree.newTerm("SO:0001183","morpholino_backbone","SO:0000348");
		tree.newTerm("SO:0001184","PNA","SO:0000348");
		tree.newTerm("SO:0001185","enzymatic","SO:0000733");
		tree.newTerm("SO:0001186","ribozymic","SO:0001185");
		tree.newTerm("SO:0001187","pseudouridylation_guide_snoRNA","SO:0000594");
		tree.newTerm("SO:0001188","LNA","SO:0000348");
		tree.newTerm("SO:0001189","LNA_oligo","SO:0001247");
		tree.newTerm("SO:0001190","TNA","SO:0000348");
		tree.newTerm("SO:0001191","TNA_oligo","SO:0001247");
		tree.newTerm("SO:0001192","GNA","SO:0000348");
		tree.newTerm("SO:0001193","GNA_oligo","SO:0001247");
		tree.newTerm("SO:0001194","R_GNA","SO:0001192");
		tree.newTerm("SO:0001195","R_GNA_oligo","SO:0001193");
		tree.newTerm("SO:0001196","S_GNA","SO:0001192");
		tree.newTerm("SO:0001197","S_GNA_oligo","SO:0001193");
		tree.newTerm("SO:0001198","ds_DNA_viral_sequence","SO:0001041");
		tree.newTerm("SO:0001199","ss_RNA_viral_sequence","SO:0001041");
		tree.newTerm("SO:0001200","negative_sense_ssRNA_viral_sequence","SO:0001199");
		tree.newTerm("SO:0001201","positive_sense_ssRNA_viral_sequence","SO:0001199");
		tree.newTerm("SO:0001202","ambisense_ssRNA_viral_sequence","SO:0001199");
		tree.newTerm("SO:0001203","RNA_polymerase_promoter","SO:0000167");
		tree.newTerm("SO:0001204","Phage_RNA_Polymerase_Promoter","SO:0001203");
		tree.newTerm("SO:0001205","SP6_RNA_Polymerase_Promoter","SO:0001204");
		tree.newTerm("SO:0001206","T3_RNA_Polymerase_Promoter","SO:0001204");
		tree.newTerm("SO:0001207","T7_RNA_Polymerase_Promoter","SO:0001204");
		tree.newTerm("SO:0001208","five_prime_EST","SO:0000345");
		tree.newTerm("SO:0001209","three_prime_EST","SO:0000345");
		tree.newTerm("SO:0001210","translational_frameshift","SO:0000836");
		tree.newTerm("SO:0001211","plus_1_translational_frameshift","SO:0001210");
		tree.newTerm("SO:0001212","plus_2_translational_frameshift","SO:0001210");
		tree.newTerm("SO:0001213","group_III_intron","SO:0000588");
		tree.newTerm("SO:0001214","noncoding_region_of_exon","SO:0000852");
		tree.newTerm("SO:0001215","coding_region_of_exon","SO:0000852");
		tree.newTerm("SO:0001216","endonuclease_spliced_intron","SO:0000188");
		tree.newTerm("SO:0001217","protein_coding_gene","SO:0000704");
		tree.newTerm("SO:0001218","transgenic_insertion","SO:0000667");
		tree.newTerm("SO:0001219","retrogene","SO:0000704");
		tree.newTerm("SO:0001220","silenced_by_RNA_interference","SO:0000893");
		tree.newTerm("SO:0001221","silenced_by_histone_modification","SO:0000893");
		tree.newTerm("SO:0001222","silenced_by_histone_methylation","SO:0001221");
		tree.newTerm("SO:0001223","silenced_by_histone_deacetylation","SO:0001221");
		tree.newTerm("SO:0001224","gene_silenced_by_RNA_interference","SO:0000127");
		tree.newTerm("SO:0001225","gene_silenced_by_histone_modification","SO:0000127");
		tree.newTerm("SO:0001226","gene_silenced_by_histone_methylation","SO:0001225");
		tree.newTerm("SO:0001227","gene_silenced_by_histone_deacetylation","SO:0001225");
		tree.newTerm("SO:0001228","dihydrouridine","SO:0001277");
		tree.newTerm("SO:0001229","pseudouridine","SO:0001277");
		tree.newTerm("SO:0001230","inosine","SO:0000250");
		tree.newTerm("SO:0001231","seven_methylguanine","SO:0000250");
		tree.newTerm("SO:0001232","ribothymidine","SO:0000250");
		tree.newTerm("SO:0001233","methylinosine","SO:0001274");
		tree.newTerm("SO:0001234","mobile","SO:0000733");
		tree.newTerm("SO:0001235","replicon","SO:0001411");
		tree.newTerm("SO:0001236","base","SO:0001411");
		tree.newTerm("SO:0001237","amino_acid","SO:0001411");
		tree.newTerm("SO:0001238","major_TSS","SO:0000315");
		tree.newTerm("SO:0001239","minor_TSS","SO:0000315");
		tree.newTerm("SO:0001240","TSS_region","SO:0000842");
		tree.newTerm("SO:0001241","encodes_alternate_transcription_start_sites","SO:0000401");
		tree.newTerm("SO:0001243","miRNA_primary_transcript_region","SO:0000835");
		tree.newTerm("SO:0001244","pre_miRNA","SO:0001243");
		tree.newTerm("SO:0001245","miRNA_stem","SO:0001243");
		tree.newTerm("SO:0001246","miRNA_loop","SO:0001243");
		tree.newTerm("SO:0001247","synthetic_oligo","SO:0000696");
		tree.newTerm("SO:0001248","assembly","SO:0001410");
		tree.newTerm("SO:0001249","fragment_assembly","SO:0001248");
		tree.newTerm("SO:0001250","fingerprint_map","SO:0001249");
		tree.newTerm("SO:0001251","STS_map","SO:0001249");
		tree.newTerm("SO:0001252","RH_map","SO:0001249");
		tree.newTerm("SO:0001253","sonicate_fragment","SO:0000143");
		tree.newTerm("SO:0001254","polyploid","SO:1000182");
		tree.newTerm("SO:0001255","autopolyploid","SO:0001254");
		tree.newTerm("SO:0001256","allopolyploid","SO:0001254");
		tree.newTerm("SO:0001257","homing_endonuclease_binding_site","SO:0000059");
		tree.newTerm("SO:0001258","octamer_motif","SO:0000713");
		tree.newTerm("SO:0001259","apicoplast_chromosome","SO:0000340");
		tree.newTerm("SO:0001260","sequence_collection");
		tree.newTerm("SO:0001261","overlapping_feature_set","SO:0000703");
		tree.newTerm("SO:0001262","overlapping_EST_set","SO:0001261");
		tree.newTerm("SO:0001263","ncRNA_gene","SO:0000704");
		tree.newTerm("SO:0001264","gRNA_gene","SO:0001263");
		tree.newTerm("SO:0001265","miRNA_gene","SO:0001263");
		tree.newTerm("SO:0001266","scRNA_gene","SO:0001263");
		tree.newTerm("SO:0001267","snoRNA_gene","SO:0001263");
		tree.newTerm("SO:0001268","snRNA_gene","SO:0001263");
		tree.newTerm("SO:0001269","SRP_RNA_gene","SO:0001263");
		tree.newTerm("SO:0001271","tmRNA_gene","SO:0001263");
		tree.newTerm("SO:0001272","tRNA_gene","SO:0001263");
		tree.newTerm("SO:0001273","modified_adenosine","SO:0000250");
		tree.newTerm("SO:0001274","modified_inosine","SO:0001230");
		tree.newTerm("SO:0001275","modified_cytidine","SO:0000250");
		tree.newTerm("SO:0001276","modified_guanosine","SO:0000250");
		tree.newTerm("SO:0001277","modified_uridine","SO:0000250");
		tree.newTerm("SO:0001278","one_methylinosine","SO:0001274");
		tree.newTerm("SO:0001279","one_two_prime_O_dimethylinosine","SO:0001274");
		tree.newTerm("SO:0001280","two_prime_O_methylinosine","SO:0001274");
		tree.newTerm("SO:0001281","three_methylcytidine","SO:0001275");
		tree.newTerm("SO:0001282","five_methylcytidine","SO:0001275");
		tree.newTerm("SO:0001283","two_prime_O_methylcytidine","SO:0001275");
		tree.newTerm("SO:0001284","two_thiocytidine","SO:0001275");
		tree.newTerm("SO:0001285","N4_acetylcytidine","SO:0001275");
		tree.newTerm("SO:0001286","five_formylcytidine","SO:0001275");
		tree.newTerm("SO:0001287","five_two_prime_O_dimethylcytidine","SO:0001275");
		tree.newTerm("SO:0001288","N4_acetyl_2_prime_O_methylcytidine","SO:0001275");
		tree.newTerm("SO:0001289","lysidine","SO:0001275");
		tree.newTerm("SO:0001290","N4_methylcytidine","SO:0001275");
		tree.newTerm("SO:0001291","N4_2_prime_O_dimethylcytidine","SO:0001275");
		tree.newTerm("SO:0001292","five_hydroxymethylcytidine","SO:0001275");
		tree.newTerm("SO:0001293","five_formyl_two_prime_O_methylcytidine","SO:0001275");
		tree.newTerm("SO:0001294","N4_N4_2_prime_O_trimethylcytidine","SO:0001275");
		tree.newTerm("SO:0001295","one_methyladenosine","SO:0001273");
		tree.newTerm("SO:0001296","two_methyladenosine","SO:0001273");
		tree.newTerm("SO:0001297","N6_methyladenosine","SO:0001273");
		tree.newTerm("SO:0001298","two_prime_O_methyladenosine","SO:0001273");
		tree.newTerm("SO:0001299","two_methylthio_N6_methyladenosine","SO:0001273");
		tree.newTerm("SO:0001300","N6_isopentenyladenosine","SO:0001273");
		tree.newTerm("SO:0001301","two_methylthio_N6_isopentenyladenosine","SO:0001273");
		tree.newTerm("SO:0001302","N6_cis_hydroxyisopentenyl_adenosine","SO:0001273");
		tree.newTerm("SO:0001303","two_methylthio_N6_cis_hydroxyisopentenyl_adenosine","SO:0001273");
		tree.newTerm("SO:0001304","N6_glycinylcarbamoyladenosine","SO:0001273");
		tree.newTerm("SO:0001305","N6_threonylcarbamoyladenosine","SO:0001273");
		tree.newTerm("SO:0001306","two_methylthio_N6_threonyl_carbamoyladenosine","SO:0001273");
		tree.newTerm("SO:0001307","N6_methyl_N6_threonylcarbamoyladenosine","SO:0001273");
		tree.newTerm("SO:0001308","N6_hydroxynorvalylcarbamoyladenosine","SO:0001273");
		tree.newTerm("SO:0001309","two_methylthio_N6_hydroxynorvalyl_carbamoyladenosine","SO:0001273");
		tree.newTerm("SO:0001310","two_prime_O_ribosyladenosine_phosphate","SO:0001273");
		tree.newTerm("SO:0001311","N6_N6_dimethyladenosine","SO:0001273");
		tree.newTerm("SO:0001312","N6_2_prime_O_dimethyladenosine","SO:0001273");
		tree.newTerm("SO:0001313","N6_N6_2_prime_O_trimethyladenosine","SO:0001273");
		tree.newTerm("SO:0001314","one_two_prime_O_dimethyladenosine","SO:0001273");
		tree.newTerm("SO:0001315","N6_acetyladenosine","SO:0001273");
		tree.newTerm("SO:0001316","seven_deazaguanosine","SO:0001276");
		tree.newTerm("SO:0001317","queuosine","SO:0001316");
		tree.newTerm("SO:0001318","epoxyqueuosine","SO:0001316");
		tree.newTerm("SO:0001319","galactosyl_queuosine","SO:0001316");
		tree.newTerm("SO:0001320","mannosyl_queuosine","SO:0001316");
		tree.newTerm("SO:0001321","seven_cyano_seven_deazaguanosine","SO:0001316");
		tree.newTerm("SO:0001322","seven_aminomethyl_seven_deazaguanosine","SO:0001316");
		tree.newTerm("SO:0001323","archaeosine","SO:0001316");
		tree.newTerm("SO:0001324","one_methylguanosine","SO:0001276");
		tree.newTerm("SO:0001325","N2_methylguanosine","SO:0001276");
		tree.newTerm("SO:0001326","seven_methylguanosine","SO:0001276");
		tree.newTerm("SO:0001327","two_prime_O_methylguanosine","SO:0001276");
		tree.newTerm("SO:0001328","N2_N2_dimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001329","N2_2_prime_O_dimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001330","N2_N2_2_prime_O_trimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001331","two_prime_O_ribosylguanosine_phosphate","SO:0001276");
		tree.newTerm("SO:0001332","wybutosine","SO:0001276");
		tree.newTerm("SO:0001333","peroxywybutosine","SO:0001276");
		tree.newTerm("SO:0001334","hydroxywybutosine","SO:0001276");
		tree.newTerm("SO:0001335","undermodified_hydroxywybutosine","SO:0001276");
		tree.newTerm("SO:0001336","wyosine","SO:0001276");
		tree.newTerm("SO:0001337","methylwyosine","SO:0001276");
		tree.newTerm("SO:0001338","N2_7_dimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001339","N2_N2_7_trimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001340","one_two_prime_O_dimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001341","four_demethylwyosine","SO:0001276");
		tree.newTerm("SO:0001342","isowyosine","SO:0001276");
		tree.newTerm("SO:0001343","N2_7_2prirme_O_trimethylguanosine","SO:0001276");
		tree.newTerm("SO:0001344","five_methyluridine","SO:0001277");
		tree.newTerm("SO:0001345","two_prime_O_methyluridine","SO:0001277");
		tree.newTerm("SO:0001346","five_two_prime_O_dimethyluridine","SO:0001277");
		tree.newTerm("SO:0001347","one_methylpseudouridine","SO:0001277");
		tree.newTerm("SO:0001348","two_prime_O_methylpseudouridine","SO:0001277");
		tree.newTerm("SO:0001349","two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001350","four_thiouridine","SO:0001277");
		tree.newTerm("SO:0001351","five_methyl_2_thiouridine","SO:0001277");
		tree.newTerm("SO:0001352","two_thio_two_prime_O_methyluridine","SO:0001277");
		tree.newTerm("SO:0001353","three_three_amino_three_carboxypropyl_uridine","SO:0001277");
		tree.newTerm("SO:0001354","five_hydroxyuridine","SO:0001277");
		tree.newTerm("SO:0001355","five_methoxyuridine","SO:0001277");
		tree.newTerm("SO:0001356","uridine_five_oxyacetic_acid","SO:0001277");
		tree.newTerm("SO:0001357","uridine_five_oxyacetic_acid_methyl_ester","SO:0001277");
		tree.newTerm("SO:0001358","five_carboxyhydroxymethyl_uridine","SO:0001277");
		tree.newTerm("SO:0001359","five_carboxyhydroxymethyl_uridine_methyl_ester","SO:0001277");
		tree.newTerm("SO:0001360","five_methoxycarbonylmethyluridine","SO:0001277");
		tree.newTerm("SO:0001361","five_methoxycarbonylmethyl_two_prime_O_methyluridine","SO:0001277");
		tree.newTerm("SO:0001362","five_methoxycarbonylmethyl_two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001363","five_aminomethyl_two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001364","five_methylaminomethyluridine","SO:0001277");
		tree.newTerm("SO:0001365","five_methylaminomethyl_two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001366","five_methylaminomethyl_two_selenouridine","SO:0001277");
		tree.newTerm("SO:0001367","five_carbamoylmethyluridine","SO:0001277");
		tree.newTerm("SO:0001368","five_carbamoylmethyl_two_prime_O_methyluridine","SO:0001277");
		tree.newTerm("SO:0001369","five_carboxymethylaminomethyluridine","SO:0001277");
		tree.newTerm("SO:0001370","five_carboxymethylaminomethyl_two_prime_O_methyluridine","SO:0001277");
		tree.newTerm("SO:0001371","five_carboxymethylaminomethyl_two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001372","three_methyluridine","SO:0001277");
		tree.newTerm("SO:0001373","one_methyl_three_three_amino_three_carboxypropyl_pseudouridine","SO:0001277");
		tree.newTerm("SO:0001374","five_carboxymethyluridine","SO:0001277");
		tree.newTerm("SO:0001375","three_two_prime_O_dimethyluridine","SO:0001277");
		tree.newTerm("SO:0001376","five_methyldihydrouridine","SO:0001277");
		tree.newTerm("SO:0001377","three_methylpseudouridine","SO:0001277");
		tree.newTerm("SO:0001378","five_taurinomethyluridine","SO:0001277");
		tree.newTerm("SO:0001379","five_taurinomethyl_two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001380","five_isopentenylaminomethyl_uridine","SO:0001277");
		tree.newTerm("SO:0001381","five_isopentenylaminomethyl_two_thiouridine","SO:0001277");
		tree.newTerm("SO:0001382","five_isopentenylaminomethyl_two_prime_O_methyluridine","SO:0001277");
		tree.newTerm("SO:0001383","histone_binding_site","SO:0001654");
		tree.newTerm("SO:0001384","CDS_fragment","SO:0000316");
		tree.newTerm("SO:0001385","modified_amino_acid_feature","SO:0001237");
		tree.newTerm("SO:0001386","modified_glycine","SO:0001385");
		tree.newTerm("SO:0001387","modified_L_alanine","SO:0001385");
		tree.newTerm("SO:0001388","modified_L_asparagine","SO:0001385");
		tree.newTerm("SO:0001389","modified_L_aspartic_acid","SO:0001385");
		tree.newTerm("SO:0001390","modified_L_cysteine","SO:0001385");
		tree.newTerm("SO:0001391","modified_L_glutamic_acid","SO:0001385");
		tree.newTerm("SO:0001392","modified_L_threonine","SO:0001385");
		tree.newTerm("SO:0001393","modified_L_tryptophan","SO:0001385");
		tree.newTerm("SO:0001394","modified_L_glutamine","SO:0001385");
		tree.newTerm("SO:0001395","modified_L_methionine","SO:0001385");
		tree.newTerm("SO:0001396","modified_L_isoleucine","SO:0001385");
		tree.newTerm("SO:0001397","modified_L_phenylalanine","SO:0001385");
		tree.newTerm("SO:0001398","modified_L_histidine","SO:0001385");
		tree.newTerm("SO:0001399","modified_L_serine","SO:0001385");
		tree.newTerm("SO:0001400","modified_L_lysine","SO:0001385");
		tree.newTerm("SO:0001401","modified_L_leucine","SO:0001385");
		tree.newTerm("SO:0001402","modified_L_selenocysteine","SO:0001385");
		tree.newTerm("SO:0001403","modified_L_valine","SO:0001385");
		tree.newTerm("SO:0001404","modified_L_proline","SO:0001385");
		tree.newTerm("SO:0001405","modified_L_tyrosine","SO:0001385");
		tree.newTerm("SO:0001406","modified_L_arginine","SO:0001385");
		tree.newTerm("SO:0001407","peptidyl","SO:0000443");
		tree.newTerm("SO:0001408","cleaved_for_gpi_anchor_region","SO:0100011");
		tree.newTerm("SO:0001409","biomaterial_region","SO:0000001");
		tree.newTerm("SO:0001410","experimental_feature","SO:0000001");
		tree.newTerm("SO:0001411","biological_region","SO:0000001");
		tree.newTerm("SO:0001412","topologically_defined_region","SO:0000001");
		tree.newTerm("SO:0001413","translocation_breakpoint","SO:0001021");
		tree.newTerm("SO:0001414","insertion_breakpoint","SO:0001021");
		tree.newTerm("SO:0001415","deletion_breakpoint","SO:0001021");
		tree.newTerm("SO:0001416","five_prime_flanking_region","SO:0000239");
		tree.newTerm("SO:0001417","three_prime_flanking_region","SO:0000239");
		tree.newTerm("SO:0001418","transcribed_fragment","SO:0001410");
		tree.newTerm("SO:0001419","cis_splice_site","SO:0000162");
		tree.newTerm("SO:0001420","trans_splice_site","SO:0000162");
		tree.newTerm("SO:0001421","splice_junction","SO:0000699");
		tree.newTerm("SO:0001422","conformational_switch","SO:0100001");
		tree.newTerm("SO:0001423","dye_terminator_read","SO:0000150");
		tree.newTerm("SO:0001424","pyrosequenced_read","SO:0000150");
		tree.newTerm("SO:0001425","ligation_based_read","SO:0000150");
		tree.newTerm("SO:0001426","polymerase_synthesis_read","SO:0000150");
		tree.newTerm("SO:0001427","cis_regulatory_frameshift_element","SO:0001679");
		tree.newTerm("SO:0001428","expressed_sequence_assembly","SO:0000353");
		tree.newTerm("SO:0001429","DNA_binding_site","SO:0001655");
		tree.newTerm("SO:0001431","cryptic_gene","SO:0000704");
		tree.newTerm("SO:0001433","three_prime_RACE_clone","SO:0000317");
		tree.newTerm("SO:0001434","cassette_pseudogene","SO:0001760");
		tree.newTerm("SO:0001435","alanine","SO:0001237");
		tree.newTerm("SO:0001436","valine","SO:0001237");
		tree.newTerm("SO:0001437","leucine","SO:0001237");
		tree.newTerm("SO:0001438","isoleucine","SO:0001237");
		tree.newTerm("SO:0001439","proline","SO:0001237");
		tree.newTerm("SO:0001440","tryptophan","SO:0001237");
		tree.newTerm("SO:0001441","phenylalanine","SO:0001237");
		tree.newTerm("SO:0001442","methionine","SO:0001237");
		tree.newTerm("SO:0001443","glycine","SO:0001237");
		tree.newTerm("SO:0001444","serine","SO:0001237");
		tree.newTerm("SO:0001445","threonine","SO:0001237");
		tree.newTerm("SO:0001446","tyrosine","SO:0001237");
		tree.newTerm("SO:0001447","cysteine","SO:0001237");
		tree.newTerm("SO:0001448","glutamine","SO:0001237");
		tree.newTerm("SO:0001449","asparagine","SO:0001237");
		tree.newTerm("SO:0001450","lysine","SO:0001237");
		tree.newTerm("SO:0001451","arginine","SO:0001237");
		tree.newTerm("SO:0001452","histidine","SO:0001237");
		tree.newTerm("SO:0001453","aspartic_acid","SO:0001237");
		tree.newTerm("SO:0001454","glutamic_acid","SO:0001237");
		tree.newTerm("SO:0001455","selenocysteine","SO:0001237");
		tree.newTerm("SO:0001456","pyrrolysine","SO:0001237");
		tree.newTerm("SO:0001457","transcribed_cluster","SO:0001410");
		tree.newTerm("SO:0001458","unigene_cluster","SO:0001457");
		tree.newTerm("SO:0001459","CRISPR","SO:0000314");
		tree.newTerm("SO:0001460","insulator_binding_site","SO:0001654");
		tree.newTerm("SO:0001461","enhancer_binding_site","SO:0001654");
		tree.newTerm("SO:0001462","contig_collection","SO:0001085","SO:0001260");
		tree.newTerm("SO:0001463","lincRNA","SO:0001877");
		tree.newTerm("SO:0001464","UST","SO:0000345");
		tree.newTerm("SO:0001465","three_prime_UST","SO:0001464");
		tree.newTerm("SO:0001466","five_prime_UST","SO:0001464");
		tree.newTerm("SO:0001467","RST","SO:0000345");
		tree.newTerm("SO:0001468","three_prime_RST","SO:0001467");
		tree.newTerm("SO:0001469","five_prime_RST","SO:0001467");
		tree.newTerm("SO:0001470","UST_match","SO:0000102");
		tree.newTerm("SO:0001471","RST_match","SO:0000102");
		tree.newTerm("SO:0001472","primer_match","SO:0000347");
		tree.newTerm("SO:0001473","miRNA_antiguide","SO:0001243");
		tree.newTerm("SO:0001474","trans_splice_junction","SO:0000699");
		tree.newTerm("SO:0001475","outron","SO:0000835");
		tree.newTerm("SO:0001476","natural_plasmid","SO:0000155","SO:0001038");
		tree.newTerm("SO:0001477","gene_trap_construct","SO:0000637");
		tree.newTerm("SO:0001478","promoter_trap_construct","SO:0000637");
		tree.newTerm("SO:0001479","enhancer_trap_construct","SO:0000637");
		tree.newTerm("SO:0001480","PAC_end","SO:0000150");
		tree.newTerm("SO:0001481","RAPD","SO:0000006");
		tree.newTerm("SO:0001482","shadow_enhancer","SO:0000165");
		tree.newTerm("SO:0001483","SNV","SO:1000002");
		tree.newTerm("SO:0001484","X_element_combinatorial_repeat","SO:0000657");
		tree.newTerm("SO:0001485","Y_prime_element","SO:0000657");
		tree.newTerm("SO:0001486","standard_draft","SO:0001499");
		tree.newTerm("SO:0001487","high_quality_draft","SO:0001499");
		tree.newTerm("SO:0001488","improved_high_quality_draft","SO:0001499");
		tree.newTerm("SO:0001489","annotation_directed_improved_draft","SO:0001499");
		tree.newTerm("SO:0001490","noncontiguous_finished","SO:0001499");
		tree.newTerm("SO:0001491","finished_genome","SO:0001499");
		tree.newTerm("SO:0001492","intronic_regulatory_region","SO:0001679");
		tree.newTerm("SO:0001493","centromere_DNA_Element_I","SO:0000330");
		tree.newTerm("SO:0001494","centromere_DNA_Element_II","SO:0000330");
		tree.newTerm("SO:0001495","centromere_DNA_Element_III","SO:0000330");
		tree.newTerm("SO:0001496","telomeric_repeat","SO:0000657");
		tree.newTerm("SO:0001497","X_element","SO:0000330");
		tree.newTerm("SO:0001498","YAC_end","SO:0000150");
		tree.newTerm("SO:0001499","whole_genome_sequence_status","SO:0000905");
		tree.newTerm("SO:0001500","heritable_phenotypic_marker","SO:0001645");
		tree.newTerm("SO:0001501","peptide_collection","SO:0001260");
		tree.newTerm("SO:0001502","high_identity_region","SO:0001410");
		tree.newTerm("SO:0001503","processed_transcript","SO:0000673");
		tree.newTerm("SO:0001504","assortment_derived_variation","SO:0000240");
		tree.newTerm("SO:0001505","reference_genome","SO:0001026");
		tree.newTerm("SO:0001506","variant_genome","SO:0001026");
		tree.newTerm("SO:0001507","variant_collection","SO:0001260");
		tree.newTerm("SO:0001508","alteration_attribute","SO:0000733");
		tree.newTerm("SO:0001509","chromosomal_variation_attribute","SO:0001508");
		tree.newTerm("SO:0001510","intrachromosomal","SO:0001509");
		tree.newTerm("SO:0001511","interchromosomal","SO:0001509");
		tree.newTerm("SO:0001512","insertion_attribute","SO:0001508");
		tree.newTerm("SO:0001513","tandem","SO:0001512");
		tree.newTerm("SO:0001514","direct","SO:0001512");
		tree.newTerm("SO:0001515","inverted","SO:0001512");
		tree.newTerm("SO:0001516","free","SO:0001523");
		tree.newTerm("SO:0001517","inversion_attribute","SO:0001508");
		tree.newTerm("SO:0001518","pericentric","SO:0001517");
		tree.newTerm("SO:0001519","paracentric","SO:0001517");
		tree.newTerm("SO:0001520","translocaton_attribute","SO:0001508");
		tree.newTerm("SO:0001521","reciprocal","SO:0001520");
		tree.newTerm("SO:0001522","insertional","SO:0001520");
		tree.newTerm("SO:0001523","duplication_attribute","SO:0001508");
		tree.newTerm("SO:0001524","chromosomally_aberrant_genome","SO:0001506");
		tree.newTerm("SO:0001525","assembly_error_correction","SO:0000413");
		tree.newTerm("SO:0001526","base_call_error_correction","SO:0000413");
		tree.newTerm("SO:0001527","peptide_localization_signal","SO:0000839");
		tree.newTerm("SO:0001528","nuclear_localization_signal","SO:0001527");
		tree.newTerm("SO:0001529","endosomal_localization_signal","SO:0001527");
		tree.newTerm("SO:0001530","lysosomal_localization_signal","SO:0001527");
		tree.newTerm("SO:0001531","nuclear_export_signal","SO:0001527");
		tree.newTerm("SO:0001532","recombination_signal_sequence","SO:0000299");
		tree.newTerm("SO:0001533","cryptic_splice_site","SO:0000162");
		tree.newTerm("SO:0001534","nuclear_rim_localization_signal","SO:0001527");
		tree.newTerm("SO:0001535","p_element","SO:0000182");
		tree.newTerm("SO:0001536","functional_variant","SO:0001060");
		tree.newTerm("SO:0001537","structural_variant","SO:0001060");
		tree.newTerm("SO:0001538","transcript_function_variant","SO:0001536");
		tree.newTerm("SO:0001539","translational_product_function_variant","SO:0001536");
		tree.newTerm("SO:0001540","level_of_transcript_variant","SO:0001538");
		tree.newTerm("SO:0001541","decreased_transcript_level_variant","SO:0001540");
		tree.newTerm("SO:0001542","increased_transcript_level_variant","SO:0001540");
		tree.newTerm("SO:0001543","transcript_processing_variant","SO:0001538");
		tree.newTerm("SO:0001544","editing_variant","SO:0001543");
		tree.newTerm("SO:0001545","polyadenylation_variant","SO:0001543");
		tree.newTerm("SO:0001546","transcript_stability_variant","SO:0001538");
		tree.newTerm("SO:0001547","decreased_transcript_stability_variant","SO:0001546");
		tree.newTerm("SO:0001548","increased_transcript_stability_variant","SO:0001546");
		tree.newTerm("SO:0001549","transcription_variant","SO:0001538");
		tree.newTerm("SO:0001550","rate_of_transcription_variant","SO:0001549");
		tree.newTerm("SO:0001551","increased_transcription_rate_variant","SO:0001550");
		tree.newTerm("SO:0001552","decreased_transcription_rate_variant","SO:0001550");
		tree.newTerm("SO:0001553","translational_product_level_variant","SO:0001539");
		tree.newTerm("SO:0001554","polypeptide_function_variant","SO:0001539");
		tree.newTerm("SO:0001555","decreased_translational_product_level","SO:0001553");
		tree.newTerm("SO:0001556","increased_translational_product_level","SO:0001553");
		tree.newTerm("SO:0001557","polypeptide_gain_of_function_variant","SO:0001554");
		tree.newTerm("SO:0001558","polypeptide_localization_variant","SO:0001554");
		tree.newTerm("SO:0001559","polypeptide_loss_of_function_variant","SO:0001554");
		tree.newTerm("SO:0001560","inactive_ligand_binding_site","SO:0001559");
		tree.newTerm("SO:0001561","polypeptide_partial_loss_of_function","SO:0001559");
		tree.newTerm("SO:0001562","polypeptide_post_translational_processing_variant","SO:0001554");
		tree.newTerm("SO:0001563","copy_number_change","SO:0001537");
		tree.newTerm("SO:0001564","gene_variant","SO:0001878");
		tree.newTerm("SO:0001565","gene_fusion","SO:0001564","SO:0001882");
		tree.newTerm("SO:0001566","regulatory_region_variant","SO:0001878");
		tree.newTerm("SO:0001567","stop_retained_variant","SO:0001590","SO:0001819");
		tree.newTerm("SO:0001568","splicing_variant","SO:0001576");
		tree.newTerm("SO:0001569","cryptic_splice_site_variant","SO:0001568");
		tree.newTerm("SO:0001570","cryptic_splice_acceptor","SO:0001569");
		tree.newTerm("SO:0001571","cryptic_splice_donor","SO:0001569");
		tree.newTerm("SO:0001572","exon_loss_variant","SO:0001568");
		tree.newTerm("SO:0001573","intron_gain_variant","SO:0001568");
		tree.newTerm("SO:0001574","splice_acceptor_variant","SO:0001629");
		tree.newTerm("SO:0001575","splice_donor_variant","SO:0001629");
		tree.newTerm("SO:0001576","transcript_variant","SO:0001564");
		tree.newTerm("SO:0001577","complex_transcript_variant","SO:0001576");
		tree.newTerm("SO:0001578","stop_lost","SO:0001590","SO:0001907","SO:0001992");
		tree.newTerm("SO:0001580","coding_sequence_variant","SO:0001791","SO:0001968");
		tree.newTerm("SO:0001582","initiator_codon_variant","SO:0001580");
		tree.newTerm("SO:0001583","missense_variant","SO:0001992");
		tree.newTerm("SO:0001585","conservative_missense_variant","SO:0001583");
		tree.newTerm("SO:0001586","non_conservative_missense_variant","SO:0001583");
		tree.newTerm("SO:0001587","stop_gained","SO:0001906","SO:0001992");
		tree.newTerm("SO:0001589","frameshift_variant","SO:0001818");
		tree.newTerm("SO:0001590","terminator_codon_variant","SO:0001580");
		tree.newTerm("SO:0001591","frame_restoring_variant","SO:0001589");
		tree.newTerm("SO:0001592","minus_1_frameshift_variant","SO:0001589");
		tree.newTerm("SO:0001593","minus_2_frameshift_variant","SO:0001589");
		tree.newTerm("SO:0001594","plus_1_frameshift_variant","SO:0001589");
		tree.newTerm("SO:0001595","plus_2_frameshift_variant","SO:0001589");
		tree.newTerm("SO:0001596","transcript_secondary_structure_variant","SO:0001576");
		tree.newTerm("SO:0001597","compensatory_transcript_secondary_structure_variant","SO:0001596");
		tree.newTerm("SO:0001598","translational_product_structure_variant","SO:0001564");
		tree.newTerm("SO:0001599","3D_polypeptide_structure_variant","SO:0001539");
		tree.newTerm("SO:0001600","complex_3D_structural_variant","SO:0001599");
		tree.newTerm("SO:0001601","conformational_change_variant","SO:0001599");
		tree.newTerm("SO:0001602","complex_change_of_translational_product_variant","SO:0001539");
		tree.newTerm("SO:0001603","polypeptide_sequence_variant","SO:0001598");
		tree.newTerm("SO:0001604","amino_acid_deletion","SO:0001603");
		tree.newTerm("SO:0001605","amino_acid_insertion","SO:0001603");
		tree.newTerm("SO:0001606","amino_acid_substitution","SO:0001603");
		tree.newTerm("SO:0001607","conservative_amino_acid_substitution","SO:0001606");
		tree.newTerm("SO:0001608","non_conservative_amino_acid_substitution","SO:0001606");
		tree.newTerm("SO:0001609","elongated_polypeptide","SO:0001603");
		tree.newTerm("SO:0001610","elongated_polypeptide_C_terminal","SO:0001609");
		tree.newTerm("SO:0001611","elongated_polypeptide_N_terminal","SO:0001609");
		tree.newTerm("SO:0001612","elongated_in_frame_polypeptide_C_terminal","SO:0001610");
		tree.newTerm("SO:0001613","elongated_out_of_frame_polypeptide_C_terminal","SO:0001610");
		tree.newTerm("SO:0001614","elongated_in_frame_polypeptide_N_terminal_elongation","SO:0001611");
		tree.newTerm("SO:0001615","elongated_out_of_frame_polypeptide_N_terminal","SO:0001611");
		tree.newTerm("SO:0001616","polypeptide_fusion","SO:0001603");
		tree.newTerm("SO:0001617","polypeptide_truncation","SO:0001603");
		tree.newTerm("SO:0001618","inactive_catalytic_site","SO:0001560");
		tree.newTerm("SO:0001619","non_coding_transcript_variant","SO:0001576");
		tree.newTerm("SO:0001620","mature_miRNA_variant","SO:0001619");
		tree.newTerm("SO:0001621","NMD_transcript_variant","SO:0001576");
		tree.newTerm("SO:0001622","UTR_variant","SO:0001791","SO:0001968");
		tree.newTerm("SO:0001623","5_prime_UTR_variant","SO:0001622");
		tree.newTerm("SO:0001624","3_prime_UTR_variant","SO:0001622");
		tree.newTerm("SO:0001626","incomplete_terminal_codon_variant","SO:0001590","SO:0001650");
		tree.newTerm("SO:0001627","intron_variant","SO:0001576");
		tree.newTerm("SO:0001628","intergenic_variant","SO:0001878");
		tree.newTerm("SO:0001629","splice_site_variant","SO:0001568","SO:0001627");
		tree.newTerm("SO:0001630","splice_region_variant","SO:0001568");
		tree.newTerm("SO:0001631","upstream_gene_variant","SO:0001628");
		tree.newTerm("SO:0001632","downstream_gene_variant","SO:0001628");
		tree.newTerm("SO:0001633","5KB_downstream_variant","SO:0001632");
		tree.newTerm("SO:0001634","500B_downstream_variant","SO:0001632");
		tree.newTerm("SO:0001635","5KB_upstream_variant","SO:0001631");
		tree.newTerm("SO:0001636","2KB_upstream_variant","SO:0001631");
		tree.newTerm("SO:0001637","rRNA_gene","SO:0001263");
		tree.newTerm("SO:0001638","piRNA_gene","SO:0001263");
		tree.newTerm("SO:0001639","RNase_P_RNA_gene","SO:0001263");
		tree.newTerm("SO:0001640","RNase_MRP_RNA_gene","SO:0001263");
		tree.newTerm("SO:0001641","lincRNA_gene","SO:0001263");
		tree.newTerm("SO:0001642","mathematically_defined_repeat","SO:0001410");
		tree.newTerm("SO:0001643","telomerase_RNA_gene","SO:0001263");
		tree.newTerm("SO:0001644","targeting_vector","SO:0000440","SO:0000804");
		tree.newTerm("SO:0001645","genetic_marker","SO:0001411");
		tree.newTerm("SO:0001646","DArT_marker","SO:0001645");
		tree.newTerm("SO:0001647","kozak_sequence","SO:0000139");
		tree.newTerm("SO:0001648","nested_transposon","SO:0000101");
		tree.newTerm("SO:0001649","nested_repeat","SO:0000657");
		tree.newTerm("SO:0001650","inframe_variant","SO:0001818");
		tree.newTerm("SO:0001653","retinoic_acid_responsive_element","SO:0000713");
		tree.newTerm("SO:0001654","nucleotide_to_protein_binding_site","SO:0000410");
		tree.newTerm("SO:0001655","nucleotide_binding_site","SO:0000409");
		tree.newTerm("SO:0001656","metal_binding_site","SO:0000409");
		tree.newTerm("SO:0001657","ligand_binding_site","SO:0000409");
		tree.newTerm("SO:0001658","nested_tandem_repeat","SO:0001649");
		tree.newTerm("SO:0001659","promoter_element","SO:0000713");
		tree.newTerm("SO:0001660","core_promoter_element","SO:0001659");
		tree.newTerm("SO:0001661","RNA_polymerase_II_TATA_box","SO:0000174");
		tree.newTerm("SO:0001662","RNA_polymerase_III_TATA_box","SO:0000174");
		tree.newTerm("SO:0001663","BREd_motif","SO:0001660");
		tree.newTerm("SO:0001664","DCE","SO:0001660");
		tree.newTerm("SO:0001665","DCE_SI","SO:0000713");
		tree.newTerm("SO:0001666","DCE_SII","SO:0000713");
		tree.newTerm("SO:0001667","DCE_SIII","SO:0000713");
		tree.newTerm("SO:0001668","proximal_promoter_element","SO:0001678");
		tree.newTerm("SO:0001669","RNApol_II_core_promoter","SO:0000170");
		tree.newTerm("SO:0001670","distal_promoter_element","SO:0001678");
		tree.newTerm("SO:0001671","bacterial_RNApol_promoter_sigma_70","SO:0000613");
		tree.newTerm("SO:0001672","bacterial_RNApol_promoter_sigma54","SO:0000613");
		tree.newTerm("SO:0001673","minus_12_signal","SO:0000713");
		tree.newTerm("SO:0001674","minus_24_signal","SO:0000713");
		tree.newTerm("SO:0001675","A_box_type_1","SO:0000619");
		tree.newTerm("SO:0001676","A_box_type_2","SO:0000619");
		tree.newTerm("SO:0001677","intermediate_element","SO:0001660");
		tree.newTerm("SO:0001678","regulatory_promoter_element","SO:0001659");
		tree.newTerm("SO:0001679","transcription_regulatory_region","SO:0005836");
		tree.newTerm("SO:0001680","translation_regulatory_region","SO:0005836");
		tree.newTerm("SO:0001681","recombination_regulatory_region","SO:0005836");
		tree.newTerm("SO:0001682","replication_regulatory_region","SO:0005836");
		tree.newTerm("SO:0001683","sequence_motif","SO:0001411");
		tree.newTerm("SO:0001684","experimental_feature_attribute","SO:0000733");
		tree.newTerm("SO:0001685","score","SO:0001684");
		tree.newTerm("SO:0001686","quality_value","SO:0001684");
		tree.newTerm("SO:0001687","restriction_enzyme_recognition_site","SO:0001954");
		tree.newTerm("SO:0001688","restriction_enzyme_cleavage_junction","SO:0000699");
		tree.newTerm("SO:0001689","five_prime_restriction_enzyme_junction","SO:0001694");
		tree.newTerm("SO:0001690","three_prime_restriction_enzyme_junction","SO:0001694");
		tree.newTerm("SO:0001691","blunt_end_restriction_enzyme_cleavage_site","SO:0001687");
		tree.newTerm("SO:0001692","sticky_end_restriction_enzyme_cleavage_site","SO:0001687");
		tree.newTerm("SO:0001693","blunt_end_restriction_enzyme_cleavage_junction","SO:0001688");
		tree.newTerm("SO:0001694","single_strand_restriction_enzyme_cleavage_site","SO:0001688");
		tree.newTerm("SO:0001695","restriction_enzyme_single_strand_overhang","SO:0001954");
		tree.newTerm("SO:0001696","experimentally_defined_binding_region","SO:0001410");
		tree.newTerm("SO:0001697","ChIP_seq_region","SO:0001696");
		tree.newTerm("SO:0001698","ASPE_primer","SO:0000112");
		tree.newTerm("SO:0001699","dCAPS_primer","SO:0000112");
		tree.newTerm("SO:0001700","histone_modification","SO:0001089","SO:0001720");
		tree.newTerm("SO:0001701","histone_methylation_site","SO:0001700");
		tree.newTerm("SO:0001702","histone_acetylation_site","SO:0001700");
		tree.newTerm("SO:0001703","H3K9_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001704","H3K14_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001705","H3K4_monomethylation_site","SO:0001734");
		tree.newTerm("SO:0001706","H3K4_trimethylation","SO:0001734");
		tree.newTerm("SO:0001707","H3K9_trimethylation_site","SO:0001736");
		tree.newTerm("SO:0001708","H3K27_monomethylation_site","SO:0001732");
		tree.newTerm("SO:0001709","H3K27_trimethylation_site","SO:0001732");
		tree.newTerm("SO:0001710","H3K79_monomethylation_site","SO:0001735");
		tree.newTerm("SO:0001711","H3K79_dimethylation_site","SO:0001735");
		tree.newTerm("SO:0001712","H3K79_trimethylation_site","SO:0001735");
		tree.newTerm("SO:0001713","H4K20_monomethylation_site","SO:0001701");
		tree.newTerm("SO:0001714","H2BK5_monomethylation_site","SO:0001701");
		tree.newTerm("SO:0001715","ISRE","SO:0001055");
		tree.newTerm("SO:0001716","histone_ubiqitination_site","SO:0001700");
		tree.newTerm("SO:0001717","H2B_ubiquitination_site","SO:0001716");
		tree.newTerm("SO:0001718","H3K18_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001719","H3K23_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001720","epigenetically_modified_region","SO:0001411");
		tree.newTerm("SO:0001722","H3K36_monomethylation_site","SO:0001733");
		tree.newTerm("SO:0001723","H3K36_dimethylation_site","SO:0001733");
		tree.newTerm("SO:0001724","H3K36_trimethylation_site","SO:0001733");
		tree.newTerm("SO:0001725","H3K4_dimethylation_site","SO:0001734");
		tree.newTerm("SO:0001726","H3K27_dimethylation_site","SO:0001732");
		tree.newTerm("SO:0001727","H3K9_monomethylation_site","SO:0001736");
		tree.newTerm("SO:0001728","H3K9_dimethylation_site","SO:0001736");
		tree.newTerm("SO:0001729","H4K16_acetylation_site","SO:0001972");
		tree.newTerm("SO:0001730","H4K5_acetylation_site","SO:0001972");
		tree.newTerm("SO:0001731","H4K8_acetylation_site","SO:0001972");
		tree.newTerm("SO:0001732","H3K27_methylation_site","SO:0001701");
		tree.newTerm("SO:0001733","H3K36_methylation_site","SO:0001701");
		tree.newTerm("SO:0001734","H3K4_methylation_site","SO:0001701");
		tree.newTerm("SO:0001735","H3K79_methylation_site","SO:0001701");
		tree.newTerm("SO:0001736","H3K9_methylation_site","SO:0001701");
		tree.newTerm("SO:0001737","histone_acylation_region","SO:0001700");
		tree.newTerm("SO:0001738","H4K_acylation_region","SO:0001737");
		tree.newTerm("SO:0001739","gene_with_non_canonical_start_codon","SO:0000704");
		tree.newTerm("SO:0001740","gene_with_start_codon_CUG","SO:0001739");
		tree.newTerm("SO:0001741","pseudogenic_gene_segment","SO:3000000");
		tree.newTerm("SO:0001742","copy_number_gain","SO:0001019");
		tree.newTerm("SO:0001743","copy_number_loss","SO:0001019");
		tree.newTerm("SO:0001744","UPD","SO:0001059");
		tree.newTerm("SO:0001745","maternal_uniparental_disomy","SO:0001744");
		tree.newTerm("SO:0001746","paternal_uniparental_disomy","SO:0001744");
		tree.newTerm("SO:0001747","open_chromatin_region","SO:0001411");
		tree.newTerm("SO:0001748","SL3_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001749","SL4_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001750","SL5_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001751","SL6_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001752","SL7_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001753","SL8_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001754","SL9_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001755","SL10_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001756","SL11_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001757","SL12_acceptor_site","SO:0000709");
		tree.newTerm("SO:0001758","duplicated_pseudogene","SO:0001760");
		tree.newTerm("SO:0001759","unitary_pseudogene","SO:0001760");
		tree.newTerm("SO:0001760","non_processed_pseudogene","SO:0000336");
		tree.newTerm("SO:0001761","variant_quality","SO:0000400");
		tree.newTerm("SO:0001762","variant_origin","SO:0001761");
		tree.newTerm("SO:0001763","variant_frequency","SO:0001761");
		tree.newTerm("SO:0001764","unique_variant","SO:0001763");
		tree.newTerm("SO:0001765","rare_variant","SO:0001763");
		tree.newTerm("SO:0001766","polymorphic_variant","SO:0001763");
		tree.newTerm("SO:0001767","common_variant","SO:0001763");
		tree.newTerm("SO:0001768","fixed_variant","SO:0001763");
		tree.newTerm("SO:0001769","variant_phenotype","SO:0001761");
		tree.newTerm("SO:0001770","benign_variant","SO:0001769");
		tree.newTerm("SO:0001771","disease_associated_variant","SO:0001769");
		tree.newTerm("SO:0001772","disease_causing_variant","SO:0001769");
		tree.newTerm("SO:0001773","lethal_variant","SO:0001536");
		tree.newTerm("SO:0001774","quantitative_variant","SO:0001769");
		tree.newTerm("SO:0001775","maternal_variant","SO:0001762");
		tree.newTerm("SO:0001776","paternal_variant","SO:0001762");
		tree.newTerm("SO:0001777","somatic_variant","SO:0001762");
		tree.newTerm("SO:0001778","germline_variant","SO:0001762");
		tree.newTerm("SO:0001779","pedigree_specific_variant","SO:0001762");
		tree.newTerm("SO:0001780","population_specific_variant","SO:0001762");
		tree.newTerm("SO:0001781","de_novo_variant","SO:0001762");
		tree.newTerm("SO:0001782","TF_binding_site_variant","SO:0001566");
		tree.newTerm("SO:0001784","complex_structural_alteration","SO:0001785","SO:1000183");
		tree.newTerm("SO:0001785","structural_alteration","SO:0001059");
		tree.newTerm("SO:0001786","loss_of_heterozygosity","SO:0001536");
		tree.newTerm("SO:0001787","splice_donor_5th_base_variant","SO:0001629");
		tree.newTerm("SO:0001788","U_box","SO:0000330");
		tree.newTerm("SO:0001789","mating_type_region","SO:0005855");
		tree.newTerm("SO:0001790","paired_end_fragment","SO:0000143");
		tree.newTerm("SO:0001791","exon_variant","SO:0001576");
		tree.newTerm("SO:0001792","non_coding_transcript_exon_variant","SO:0001619","SO:0001791");
		tree.newTerm("SO:0001793","clone_end","SO:0000150");
		tree.newTerm("SO:0001794","point_centromere","SO:0000577");
		tree.newTerm("SO:0001795","regional_centromere","SO:0000577");
		tree.newTerm("SO:0001796","regional_centromere_central_core","SO:0000330");
		tree.newTerm("SO:0001797","centromeric_repeat","SO:0000657");
		tree.newTerm("SO:0001798","regional_centromere_inner_repeat_region","SO:0001797");
		tree.newTerm("SO:0001799","regional_centromere_outer_repeat_region","SO:0001797");
		tree.newTerm("SO:0001800","tasiRNA","SO:0000655");
		tree.newTerm("SO:0001801","tasiRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0001802","increased_polyadenylation_variant","SO:0001545");
		tree.newTerm("SO:0001803","decreased_polyadenylation_variant","SO:0001545");
		tree.newTerm("SO:0001804","DDB_box","SO:0001093");
		tree.newTerm("SO:0001805","destruction_box","SO:0100017");
		tree.newTerm("SO:0001806","ER_retention_signal","SO:0001527");
		tree.newTerm("SO:0001807","KEN_box","SO:0100017");
		tree.newTerm("SO:0001808","mitochondrial_targeting_signal","SO:0001527");
		tree.newTerm("SO:0001809","signal_anchor","SO:0000418");
		tree.newTerm("SO:0001810","PIP_box","SO:0001093");
		tree.newTerm("SO:0001811","phosphorylation_site","SO:0001089");
		tree.newTerm("SO:0001812","transmembrane_helix","SO:0001114");
		tree.newTerm("SO:0001813","vacuolar_sorting_signal","SO:0001527");
		tree.newTerm("SO:0001814","coding_variant_quality","SO:0001761");
		tree.newTerm("SO:0001815","synonymous","SO:0001814");
		tree.newTerm("SO:0001816","non_synonymous","SO:0001814");
		tree.newTerm("SO:0001817","inframe","SO:0000863");
		tree.newTerm("SO:0001818","protein_altering_variant","SO:0001580");
		tree.newTerm("SO:0001819","synonymous_variant","SO:0001580");
		tree.newTerm("SO:0001820","inframe_indel","SO:0001650");
		tree.newTerm("SO:0001821","inframe_insertion","SO:0001820","SO:0001908");
		tree.newTerm("SO:0001822","inframe_deletion","SO:0001820","SO:0001906");
		tree.newTerm("SO:0001823","conservative_inframe_insertion","SO:0001821");
		tree.newTerm("SO:0001824","disruptive_inframe_insertion","SO:0001821");
		tree.newTerm("SO:0001825","conservative_inframe_deletion","SO:0001822");
		tree.newTerm("SO:0001826","disruptive_inframe_deletion","SO:0001822");
		tree.newTerm("SO:0001827","mRNA_read","SO:0000150");
		tree.newTerm("SO:0001828","genomic_DNA_read","SO:0000150");
		tree.newTerm("SO:0001829","mRNA_contig","SO:0000149");
		tree.newTerm("SO:0001830","AFLP_fragment","SO:0000006");
		tree.newTerm("SO:0001831","protein_hmm_match","SO:0000349");
		tree.newTerm("SO:0001832","immunoglobulin_region","SO:0000839");
		tree.newTerm("SO:0001833","V_region","SO:0001832");
		tree.newTerm("SO:0001834","C_region","SO:0001832");
		tree.newTerm("SO:0001835","N_region","SO:0000301");
		tree.newTerm("SO:0001836","S_region","SO:0000301");
		tree.newTerm("SO:0001837","mobile_element_insertion","SO:0000667");
		tree.newTerm("SO:0001838","novel_sequence_insertion","SO:0000667");
		tree.newTerm("SO:0001839","CSL_response_element","SO:0001659");
		tree.newTerm("SO:0001840","GATA_box","SO:0001660");
		tree.newTerm("SO:0001841","polymorphic_pseudogene","SO:0000336");
		tree.newTerm("SO:0001842","AP_1_binding_site","SO:0001659");
		tree.newTerm("SO:0001843","CRE","SO:0001659");
		tree.newTerm("SO:0001844","CuRE","SO:0001659");
		tree.newTerm("SO:0001845","DRE","SO:0001659");
		tree.newTerm("SO:0001846","FLEX_element","SO:0001659");
		tree.newTerm("SO:0001847","forkhead_motif","SO:0001659");
		tree.newTerm("SO:0001848","homol_D_box","SO:0001660");
		tree.newTerm("SO:0001849","homol_E_box","SO:0001659");
		tree.newTerm("SO:0001850","HSE","SO:0001659");
		tree.newTerm("SO:0001851","iron_repressed_GATA_element","SO:0001840");
		tree.newTerm("SO:0001852","mating_type_M_box","SO:0001659");
		tree.newTerm("SO:0001853","androgen_response_element","SO:0000713");
		tree.newTerm("SO:0001854","smFISH_probe","SO:0000051");
		tree.newTerm("SO:0001855","MCB","SO:0001659");
		tree.newTerm("SO:0001856","CCAAT_motif","SO:0001659");
		tree.newTerm("SO:0001857","Ace2_UAS","SO:0001659");
		tree.newTerm("SO:0001858","TR_box","SO:0001659");
		tree.newTerm("SO:0001859","STREP_motif","SO:0001659");
		tree.newTerm("SO:0001860","rDNA_intergenic_spacer_element","SO:0000713");
		tree.newTerm("SO:0001861","sterol_regulatory_element","SO:0001659");
		tree.newTerm("SO:0001862","GT_dinucleotide_repeat","SO:0000290");
		tree.newTerm("SO:0001863","GTT_trinucleotide_repeat","SO:0000291");
		tree.newTerm("SO:0001864","Sap1_recognition_motif","SO:0000713");
		tree.newTerm("SO:0001865","CDRE_motif","SO:0001659");
		tree.newTerm("SO:0001866","BAC_read_contig","SO:0000149");
		tree.newTerm("SO:0001867","candidate_gene","SO:0000704");
		tree.newTerm("SO:0001868","positional_candidate_gene","SO:0001867");
		tree.newTerm("SO:0001869","functional_candidate_gene","SO:0001867");
		tree.newTerm("SO:0001870","enhancerRNA","SO:0000655");
		tree.newTerm("SO:0001871","PCB","SO:0001659");
		tree.newTerm("SO:0001872","rearrangement_region","SO:0001411","SO:0001785");
		tree.newTerm("SO:0001873","interchromosomal_breakpoint","SO:0001021");
		tree.newTerm("SO:0001874","intrachromosomal_breakpoint","SO:0001021");
		tree.newTerm("SO:0001875","unassigned_supercontig","SO:0000148");
		tree.newTerm("SO:0001876","partial_genomic_sequence_assembly","SO:0000353");
		tree.newTerm("SO:0001877","lnc_RNA","SO:0000655");
		tree.newTerm("SO:0001878","feature_variant","SO:0001537");
		tree.newTerm("SO:0001879","feature_ablation","SO:0001537");
		tree.newTerm("SO:0001880","feature_amplification","SO:0001537");
		tree.newTerm("SO:0001881","feature_translocation","SO:0001537");
		tree.newTerm("SO:0001882","feature_fusion","SO:0001537");
		tree.newTerm("SO:0001883","transcript_translocation","SO:0001881");
		tree.newTerm("SO:0001884","regulatory_region_translocation","SO:0001881");
		tree.newTerm("SO:0001885","TFBS_translocation","SO:0001884");
		tree.newTerm("SO:0001886","transcript_fusion","SO:0001882");
		tree.newTerm("SO:0001887","regulatory_region_fusion","SO:0001882");
		tree.newTerm("SO:0001888","TFBS_fusion","SO:0001887");
		tree.newTerm("SO:0001889","transcript_amplification","SO:0001880");
		tree.newTerm("SO:0001890","transcript_regulatory_region_fusion","SO:0001882");
		tree.newTerm("SO:0001891","regulatory_region_amplification","SO:0001880");
		tree.newTerm("SO:0001892","TFBS_amplification","SO:0001891");
		tree.newTerm("SO:0001893","transcript_ablation","SO:0001879");
		tree.newTerm("SO:0001894","regulatory_region_ablation","SO:0001879");
		tree.newTerm("SO:0001895","TFBS_ablation","SO:0001894");
		tree.newTerm("SO:0001896","transposable_element_CDS","SO:0000316");
		tree.newTerm("SO:0001897","transposable_element_pseudogene","SO:0000336");
		tree.newTerm("SO:0001898","dg_repeat","SO:0001797");
		tree.newTerm("SO:0001899","dh_repeat","SO:0001797");
		tree.newTerm("SO:0001900","M26_binding_site","SO:0000713");
		tree.newTerm("SO:0001901","AACCCT_box","SO:0001660");
		tree.newTerm("SO:0001902","splice_region","SO:0000835");
		tree.newTerm("SO:0001903","intronic_lncRNA","SO:0001877");
		tree.newTerm("SO:0001904","antisense_lncRNA","SO:0001877");
		tree.newTerm("SO:0001905","regional_centromere_outer_repeat_transcript","SO:0000185");
		tree.newTerm("SO:0001906","feature_truncation","SO:0001878");
		tree.newTerm("SO:0001907","feature_elongation","SO:0001878");
		tree.newTerm("SO:0001908","internal_feature_elongation","SO:0001907");
		tree.newTerm("SO:0001909","frameshift_elongation","SO:0001589","SO:0001908");
		tree.newTerm("SO:0001910","frameshift_truncation","SO:0001589","SO:0001906");
		tree.newTerm("SO:0001911","copy_number_increase","SO:0001563");
		tree.newTerm("SO:0001912","copy_number_decrease","SO:0001563");
		tree.newTerm("SO:0001913","bacterial_RNApol_promoter_sigma_ecf","SO:0000613");
		tree.newTerm("SO:0001914","rDNA_replication_fork_barrier","SO:0000713");
		tree.newTerm("SO:0001915","transcription_start_cluster","SO:0001410");
		tree.newTerm("SO:0001916","CAGE_tag","SO:0000324");
		tree.newTerm("SO:0001917","CAGE_cluster","SO:0001915");
		tree.newTerm("SO:0001918","5_methylcytosine","SO:0000114");
		tree.newTerm("SO:0001919","4_methylcytosine","SO:0000114");
		tree.newTerm("SO:0001920","N6_methyladenine","SO:0000161");
		tree.newTerm("SO:0001921","mitochondrial_contig","SO:0000149");
		tree.newTerm("SO:0001922","mitochondrial_supercontig","SO:0000148");
		tree.newTerm("SO:0001923","TERRA","SO:0001927");
		tree.newTerm("SO:0001924","ARRET","SO:0001927");
		tree.newTerm("SO:0001925","ARIA","SO:0001927");
		tree.newTerm("SO:0001926","anti_ARRET","SO:0001927");
		tree.newTerm("SO:0001927","telomeric_transcript","SO:0000655");
		tree.newTerm("SO:0001928","distal_duplication","SO:1000035");
		tree.newTerm("SO:0001929","mitochondrial_DNA_read","SO:0000150");
		tree.newTerm("SO:0001930","chloroplast_DNA_read","SO:0000150");
		tree.newTerm("SO:0001931","consensus_gDNA","SO:0000994");
		tree.newTerm("SO:0001932","restriction_enzyme_five_prime_single_strand_overhang","SO:0001695");
		tree.newTerm("SO:0001933","restriction_enzyme_three_prime_single_strand_overhang","SO:0001695");
		tree.newTerm("SO:0001934","monomeric_repeat","SO:0000705");
		tree.newTerm("SO:0001935","H3K20_trimethylation_site","SO:0001701");
		tree.newTerm("SO:0001936","H3K36_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001937","H2BK12_acetylation_site","SO:0002143");
		tree.newTerm("SO:0001938","H2AK5_acetylation_site","SO:0002142");
		tree.newTerm("SO:0001939","H4K12_acetylation_site","SO:0001972");
		tree.newTerm("SO:0001940","H2BK120_acetylation_site","SO:0002143");
		tree.newTerm("SO:0001941","H4K91_acetylation_site","SO:0001972");
		tree.newTerm("SO:0001942","H2BK20_acetylation_site","SO:0002143");
		tree.newTerm("SO:0001943","H3K4_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001944","H2AK9_acetylation_site","SO:0002142");
		tree.newTerm("SO:0001945","H3K56_acetylation_site","SO:0001973");
		tree.newTerm("SO:0001946","H2BK15_acetylation_site","SO:0002143");
		tree.newTerm("SO:0001947","H3R2_monomethylation_site","SO:0001701");
		tree.newTerm("SO:0001948","H3R2_dimethylation_site","SO:0001701");
		tree.newTerm("SO:0001949","H4R3_dimethylation_site","SO:0001701");
		tree.newTerm("SO:0001950","H4K4_trimethylation_site","SO:0001701");
		tree.newTerm("SO:0001951","H3K23_dimethylation_site","SO:0001701");
		tree.newTerm("SO:0001952","promoter_flanking_region","SO:0001055");
		tree.newTerm("SO:0001953","restriction_enzyme_assembly_scar","SO:0001954");
		tree.newTerm("SO:0001954","restriction_enzyme_region","SO:0001411");
		tree.newTerm("SO:0001955","protein_stability_element","SO:0000839");
		tree.newTerm("SO:0001956","protease_site","SO:0000839");
		tree.newTerm("SO:0001958","lariat_intron","SO:0000188");
		tree.newTerm("SO:0001959","TCT_motif","SO:0001660");
		tree.newTerm("SO:0001960","5_hydroxymethylcytosine","SO:0000114");
		tree.newTerm("SO:0001961","5_formylcytosine","SO:0001963");
		tree.newTerm("SO:0001962","modified_adenine","SO:0000305");
		tree.newTerm("SO:0001963","modified_cytosine","SO:0000305");
		tree.newTerm("SO:0001964","modified_guanine","SO:0000305");
		tree.newTerm("SO:0001965","8_oxoguanine","SO:0001964");
		tree.newTerm("SO:0001966","5_carboxylcytosine","SO:0001963");
		tree.newTerm("SO:0001967","8_oxoadenine","SO:0001962");
		tree.newTerm("SO:0001968","coding_transcript_variant","SO:0001576");
		tree.newTerm("SO:0001969","coding_transcript_intron_variant","SO:0001627","SO:0001968");
		tree.newTerm("SO:0001970","non_coding_transcript_intron_variant","SO:0001619","SO:0001627");
		tree.newTerm("SO:0001971","zinc_finger_binding_site","SO:0001429");
		tree.newTerm("SO:0001972","histone_4_acetylation_site","SO:0001702");
		tree.newTerm("SO:0001973","histone_3_acetylation_site","SO:0001702");
		tree.newTerm("SO:0001974","CTCF_binding_site","SO:0001659");
		tree.newTerm("SO:0001975","five_prime_sticky_end_restriction_enzyme_cleavage_site","SO:0001692");
		tree.newTerm("SO:0001976","three_prime_sticky_end_restriction_enzyme_cleavage_site","SO:0001692");
		tree.newTerm("SO:0001977","ribonuclease_site","SO:0000833");
		tree.newTerm("SO:0001978","signature","SO:0000804");
		tree.newTerm("SO:0001979","RNA_stability_element","SO:0000715");
		tree.newTerm("SO:0001980","G_box","SO:0001678");
		tree.newTerm("SO:0001981","L_box","SO:0001678");
		tree.newTerm("SO:0001982","I-box","SO:0001678");
		tree.newTerm("SO:0001983","5_prime_UTR_premature_start_codon_variant","SO:0001623");
		tree.newTerm("SO:0001984","silent_mating_type_cassette_array","SO:0005854");
		tree.newTerm("SO:0001985","Okazaki_fragment","SO:0001411");
		tree.newTerm("SO:0001986","upstream_transcript_variant","SO:0001628");
		tree.newTerm("SO:0001987","downstream_transcript_variant","SO:0001628");
		tree.newTerm("SO:0001988","5_prime_UTR_premature_start_codon_gain_variant","SO:0001983");
		tree.newTerm("SO:0001989","5_prime_UTR_premature_start_codon_loss_variant","SO:0001983");
		tree.newTerm("SO:0001990","five_prime_UTR_premature_start_codon_location_variant","SO:0001983");
		tree.newTerm("SO:0001991","consensus_AFLP_fragment","SO:0000994");
		tree.newTerm("SO:0001992","nonsynonymous_variant","SO:0001650");
		tree.newTerm("SO:0001993","extended_cis_splice_site","SO:0001419");
		tree.newTerm("SO:0001994","intron_base_5","SO:0001014");
		tree.newTerm("SO:0001995","extended_intronic_splice_region_variant","SO:0001568");
		tree.newTerm("SO:0001996","extended_intronic_splice_region","SO:0001014");
		tree.newTerm("SO:0001997","subtelomere","SO:0000628");
		tree.newTerm("SO:0001998","sgRNA","SO:0000696");
		tree.newTerm("SO:0001999","mating_type_region_motif","SO:0000713");
		tree.newTerm("SO:0002001","Y_region","SO:0001999");
		tree.newTerm("SO:0002002","Z1_region","SO:0001999");
		tree.newTerm("SO:0002003","Z2_region","SO:0001999");
		tree.newTerm("SO:0002004","ARS_consensus_sequence","SO:0000713");
		tree.newTerm("SO:0002005","DSR_motif","SO:0000713");
		tree.newTerm("SO:0002006","zinc_repressed_element","SO:0001659");
		tree.newTerm("SO:0002007","MNV","SO:1000002");
		tree.newTerm("SO:0002008","rare_amino_acid_variant","SO:0001586");
		tree.newTerm("SO:0002009","selenocysteine_loss","SO:0002008");
		tree.newTerm("SO:0002010","pyrrolysine_loss","SO:0002008");
		tree.newTerm("SO:0002011","intragenic_variant","SO:0001576");
		tree.newTerm("SO:0002012","start_lost","SO:0001582","SO:0001992");
		tree.newTerm("SO:0002013","5_prime_UTR_truncation","SO:0001623");
		tree.newTerm("SO:0002014","5_prime_UTR_elongation","SO:0001623");
		tree.newTerm("SO:0002015","3_prime_UTR_truncation","SO:0001624");
		tree.newTerm("SO:0002016","3_prime_UTR_elongation","SO:0001624");
		tree.newTerm("SO:0002017","conserved_intergenic_variant","SO:0001628");
		tree.newTerm("SO:0002018","conserved_intron_variant","SO:0001627");
		tree.newTerm("SO:0002019","start_retained_variant","SO:0001582","SO:0001819");
		tree.newTerm("SO:0002020","boundary_element","SO:0000713");
		tree.newTerm("SO:0002021","mating_type_region_replication_fork_barrier","SO:0000713");
		tree.newTerm("SO:0002022","priRNA","SO:0000655");
		tree.newTerm("SO:0002023","multiplexing_sequence_identifier","SO:0000324");
		tree.newTerm("SO:0002024","W_region","SO:0001999");
		tree.newTerm("SO:0002025","cis_acting_homologous_chromosome_pairing_region","SO:0000713");
		tree.newTerm("SO:0002026","intein_encoding_region","SO:0000842");
		tree.newTerm("SO:0002027","uORF","SO:0000236");
		tree.newTerm("SO:0002028","sORF","SO:0000236");
		tree.newTerm("SO:0002029","tnaORF","SO:0000236");
		tree.newTerm("SO:0002030","X_region","SO:0001999");
		tree.newTerm("SO:0002031","shRNA","SO:0000655");
		tree.newTerm("SO:0002032","moR","SO:0000370");
		tree.newTerm("SO:0002033","loR","SO:0000370");
		tree.newTerm("SO:0002034","miR_encoding_snoRNA_primary_transcript","SO:0000232");
		tree.newTerm("SO:0002035","lncRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0002036","miR_encoding_lncRNA_primary_transcript","SO:0002035");
		tree.newTerm("SO:0002037","miR_encoding_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0002038","shRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0002039","miR_encoding_shRNA_primary_transcript","SO:0002038");
		tree.newTerm("SO:0002040","vaultRNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0002041","miR_encoding_vaultRNA_primary_transcript","SO:0002040");
		tree.newTerm("SO:0002042","Y_RNA_primary_transcript","SO:0000483");
		tree.newTerm("SO:0002043","miR_encoding_Y_RNA_primary_transcript","SO:0002042");
		tree.newTerm("SO:0002044","TCS_element","SO:0001659");
		tree.newTerm("SO:0002045","pheromone_response_element","SO:0000235","SO:0001659");
		tree.newTerm("SO:0002046","FRE","SO:0001461");
		tree.newTerm("SO:0002047","transcription_pause_site","SO:0001679");
		tree.newTerm("SO:0002048","disabled_reading_frame","SO:0000717");
		tree.newTerm("SO:0002049","H3K27_acetylation_site","SO:0001973");
		tree.newTerm("SO:0002050","constitutive_promoter","SO:0000167");
		tree.newTerm("SO:0002051","inducible_promoter","SO:0000167");
		tree.newTerm("SO:0002052","dominant_negative_variant","SO:0001536");
		tree.newTerm("SO:0002053","gain_of_function_variant","SO:0001536");
		tree.newTerm("SO:0002054","loss_of_function_variant","SO:0001536");
		tree.newTerm("SO:0002055","null_mutation","SO:0001536");
		tree.newTerm("SO:0002056","intronic_splicing_silencer","SO:0000625");
		tree.newTerm("SO:0002058","exonic_splicing_silencer","SO:0000625","SO:0000852");
		tree.newTerm("SO:0002059","recombination_enhancer","SO:0001681");
		tree.newTerm("SO:0002060","interchromosomal_translocation","SO:1000044");
		tree.newTerm("SO:0002061","intrachromosomal_translocation","SO:1000044");
		tree.newTerm("SO:0002062","complex_chromosomal_rearrangement","SO:0001784");
		tree.newTerm("SO:0002063","Alu_insertion","SO:0001837");
		tree.newTerm("SO:0002064","LINE1_insertion","SO:0001837");
		tree.newTerm("SO:0002065","SVA_insertion","SO:0001837");
		tree.newTerm("SO:0002066","mobile_element_deletion","SO:0000159");
		tree.newTerm("SO:0002067","HERV_deletion","SO:0002066");
		tree.newTerm("SO:0002068","SVA_deletion","SO:0002066");
		tree.newTerm("SO:0002069","LINE1_deletion","SO:0002066");
		tree.newTerm("SO:0002070","Alu_deletion","SO:0002066");
		tree.newTerm("SO:0002071","CDS_supported_by_peptide_spectrum_match","SO:1001251");
		tree.newTerm("SO:0002072","sequence_comparison","SO:0000110");
		tree.newTerm("SO:0002073","no_sequence_alteration","SO:0002072");
		tree.newTerm("SO:0002074","intergenic_1kb_variant","SO:0001628");
		tree.newTerm("SO:0002075","incomplete_transcript_variant","SO:0001576");
		tree.newTerm("SO:0002076","incomplete_transcript_3UTR_variant","SO:0002075");
		tree.newTerm("SO:0002077","incomplete_transcript_5UTR_variant","SO:0002075");
		tree.newTerm("SO:0002078","incomplete_transcript_intronic_variant","SO:0002075");
		tree.newTerm("SO:0002079","incomplete_transcript_splice_region_variant","SO:0002075");
		tree.newTerm("SO:0002080","incomplete_transcript_exonic_variant","SO:0002075");
		tree.newTerm("SO:0002081","incomplete_transcript_CDS","SO:0002075");
		tree.newTerm("SO:0002082","incomplete_transcript_coding_splice_variant","SO:0002079");
		tree.newTerm("SO:0002083","2KB_downstream_variant","SO:0001632");
		tree.newTerm("SO:0002084","exonic_splice_region_variant","SO:0001630");
		tree.newTerm("SO:0002085","unidirectional_gene_fusion","SO:0001565");
		tree.newTerm("SO:0002086","bidirectional_gene_fusion","SO:0001565");
		tree.newTerm("SO:0002087","pseudogenic_CDS","SO:0000462");
		tree.newTerm("SO:0002088","non_coding_transcript_splice_region_variant","SO:0001619","SO:0001630");
		tree.newTerm("SO:0002089","3_prime_UTR_exon_variant","SO:0001624");
		tree.newTerm("SO:0002090","3_prime_UTR_intron_variant","SO:0001624","SO:0001969");
		tree.newTerm("SO:0002091","5_prime_UTR_intron_variant","SO:0001623","SO:0001969");
		tree.newTerm("SO:0002092","5_prime_UTR_exon_variant","SO:0001623");
		tree.newTerm("SO:0002093","structural_interaction_variant","SO:0001599");
		tree.newTerm("SO:0002094","non_allelic_homologous_recombination_region","SO:0000339");
		tree.newTerm("SO:0002095","scaRNA","SO:0000655");
		tree.newTerm("SO:0002096","short_tandem_repeat_variation","SO:0000248");
		tree.newTerm("SO:0002097","vertebrate_immune_system_pseudogene","SO:0000336");
		tree.newTerm("SO:0002098","immunoglobulin_pseudogene","SO:0002097");
		tree.newTerm("SO:0002099","T_cell_receptor_pseudogene","SO:0002097");
		tree.newTerm("SO:0002100","IG_C_pseudogene","SO:0002098");
		tree.newTerm("SO:0002101","IG_J_pseudogene","SO:0002098");
		tree.newTerm("SO:0002102","IG_V_pseudogene","SO:0002098");
		tree.newTerm("SO:0002103","TR_V_pseudogene","SO:0002099");
		tree.newTerm("SO:0002104","TR_J_pseudogene","SO:0002099");
		tree.newTerm("SO:0002105","translated_processed_pseudogene","SO:0000043");
		tree.newTerm("SO:0002106","translated_unprocessed_pseudogene","SO:0001760");
		tree.newTerm("SO:0002107","transcribed_unprocessed_pseudogene","SO:0001760");
		tree.newTerm("SO:0002108","transcribed_unitary_pseudogene","SO:0001759");
		tree.newTerm("SO:0002109","transcribed_processed_pseudogene","SO:0000043");
		tree.newTerm("SO:0002110","polymorphic_pseudogene_with_retained_intron","SO:0001841");
		tree.newTerm("SO:0002111","pseudogene_processed_transcript","SO:0001503");
		tree.newTerm("SO:0002112","coding_transcript_with_retained_intron","SO:0000120");
		tree.newTerm("SO:0002113","lncRNA_with_retained_intron","SO:0002035");
		tree.newTerm("SO:0002114","NMD_transcript","SO:0000120");
		tree.newTerm("SO:0002115","pseudogenic_transcript_with_retained_intron","SO:0000185");
		tree.newTerm("SO:0002116","polymorphic_pseudogene_processed_transcript","SO:0002111");
		tree.newTerm("SO:0002118","NMD_polymorphic_pseudogene_transcript","SO:0002114");
		tree.newTerm("SO:0002119","allelic_frequency","SO:0001763");
		tree.newTerm("SO:0002120","three_prime_overlapping_ncrna","SO:0000655");
		tree.newTerm("SO:0002121","vertebrate_immune_system_gene","SO:0000704");
		tree.newTerm("SO:0002122","immunoglobulin_gene","SO:0002121");
		tree.newTerm("SO:0002123","IG_C_gene","SO:0002122");
		tree.newTerm("SO:0002124","IG_D_gene","SO:0002122");
		tree.newTerm("SO:0002125","IG_J_gene","SO:0002122");
		tree.newTerm("SO:0002126","IG_V_gene","SO:0002122");
		tree.newTerm("SO:0002127","lncRNA_gene","SO:0001263");
		tree.newTerm("SO:0002128","mt_rRNA","SO:0000252");
		tree.newTerm("SO:0002129","mt_tRNA","SO:0000253");
		tree.newTerm("SO:0002130","NSD_transcript","SO:0000234");
		tree.newTerm("SO:0002131","sense_intronic_ncRNA","SO:0000655");
		tree.newTerm("SO:0002132","sense_overlap_ncRNA","SO:0000655");
		tree.newTerm("SO:0002133","T_cell_receptor_gene","SO:0002121");
		tree.newTerm("SO:0002134","TR_C_Gene","SO:0002133");
		tree.newTerm("SO:0002135","TR_D_Gene","SO:0002133");
		tree.newTerm("SO:0002136","TR_J_Gene","SO:0002133");
		tree.newTerm("SO:0002137","TR_V_Gene","SO:0002133");
		tree.newTerm("SO:0002138","predicted_transcript","SO:0000673");
		tree.newTerm("SO:0002139","unconfirmed_transcript","SO:0002138");
		tree.newTerm("SO:0002140","early_origin_of_replication","SO:0000296");
		tree.newTerm("SO:0002141","late_origin_of_replication","SO:0000296");
		tree.newTerm("SO:0002142","histone_2A_acetylation_site","SO:0001702");
		tree.newTerm("SO:0002143","histone_2B_acetylation_site","SO:0001702");
		tree.newTerm("SO:0002144","histone_2AZ_acetylation_site","SO:0002142");
		tree.newTerm("SO:0002145","H2AZK4_acetylation_site","SO:0002144");
		tree.newTerm("SO:0002146","H2AZK7_acetylation_site","SO:0002144");
		tree.newTerm("SO:0002147","H2AZK11_acetylation_site","SO:0002144");
		tree.newTerm("SO:0002148","H2AZK13_acetylation_site","SO:0002144");
		tree.newTerm("SO:0002149","H2AZK15_acetylation_site","SO:0002144");
		tree.newTerm("SO:0002150","AUG_initiated_uORF","SO:0002027");
		tree.newTerm("SO:0002151","non_AUG_initiated_uORF","SO:0002027");
		tree.newTerm("SO:0002152","genic_downstream_transcript_variant","SO:0001564");
		tree.newTerm("SO:0002153","genic_upstream_transcript_variant","SO:0001564");
		tree.newTerm("SO:0002154","mitotic_recombination_region","SO:0000298");
		tree.newTerm("SO:0002155","meiotic_recombination_region","SO:0000298");
		tree.newTerm("SO:0002156","CArG_box","SO:0001659");
		tree.newTerm("SO:0002157","Mat2P","SO:0001984");
		tree.newTerm("SO:0002158","Mat3M","SO:0001984");
		tree.newTerm("SO:0002159","SHP_box","SO:0100017");
		tree.newTerm("SO:0005836","regulatory_region","SO:0000831");
		tree.newTerm("SO:0005837","U14_snoRNA_primary_transcript","SO:0000232");
		tree.newTerm("SO:0005841","methylation_guide_snoRNA","SO:0000593");
		tree.newTerm("SO:0005843","rRNA_cleavage_RNA","SO:0000655");
		tree.newTerm("SO:0005845","exon_of_single_exon_gene","SO:0000147");
		tree.newTerm("SO:0005847","cassette_array_member","SO:0005848");
		tree.newTerm("SO:0005848","gene_cassette_member","SO:0000081");
		tree.newTerm("SO:0005849","gene_subarray_member","SO:0000081");
		tree.newTerm("SO:0005850","primer_binding_site","SO:0001655");
		tree.newTerm("SO:0005851","gene_array","SO:0005855");
		tree.newTerm("SO:0005852","gene_subarray","SO:0005855");
		tree.newTerm("SO:0005853","gene_cassette","SO:0000704");
		tree.newTerm("SO:0005854","gene_cassette_array","SO:0005855");
		tree.newTerm("SO:0005855","gene_group","SO:0001411");
		tree.newTerm("SO:0005856","selenocysteine_tRNA_primary_transcript","SO:0000210");
		tree.newTerm("SO:0005857","selenocysteinyl_tRNA","SO:0000253");
		tree.newTerm("SO:0005858","syntenic_region","SO:0000330");
		tree.newTerm("SO:0100001","biochemical_region_of_peptide","SO:0001067");
		tree.newTerm("SO:0100002","molecular_contact_region","SO:0100001");
		tree.newTerm("SO:0100003","intrinsically_unstructured_polypeptide_region","SO:0001070");
		tree.newTerm("SO:0100004","catmat_left_handed_three","SO:0001078");
		tree.newTerm("SO:0100005","catmat_left_handed_four","SO:0001078");
		tree.newTerm("SO:0100006","catmat_right_handed_three","SO:0001078");
		tree.newTerm("SO:0100007","catmat_right_handed_four","SO:0001078");
		tree.newTerm("SO:0100008","alpha_beta_motif","SO:0001078");
		tree.newTerm("SO:0100009","lipoprotein_signal_peptide","SO:0100011");
		tree.newTerm("SO:0100010","no_output","SO:0000703");
		tree.newTerm("SO:0100011","cleaved_peptide_region","SO:0000839");
		tree.newTerm("SO:0100012","peptide_coil","SO:0001078");
		tree.newTerm("SO:0100013","hydrophobic_region_of_peptide","SO:0000839");
		tree.newTerm("SO:0100014","n_terminal_region","SO:0100011");
		tree.newTerm("SO:0100015","c_terminal_region","SO:0100011");
		tree.newTerm("SO:0100016","central_hydrophobic_region_of_signal_peptide","SO:0100011");
		tree.newTerm("SO:0100017","polypeptide_conserved_motif","SO:0001067");
		tree.newTerm("SO:0100018","polypeptide_binding_motif","SO:0100001");
		tree.newTerm("SO:0100019","polypeptide_catalytic_motif","SO:0100001");
		tree.newTerm("SO:0100020","polypeptide_DNA_contact","SO:0001429","SO:0100002");
		tree.newTerm("SO:0100021","polypeptide_conserved_region","SO:0000839");
		tree.newTerm("SO:1000002","substitution","SO:0001059","SO:0001411");
		tree.newTerm("SO:1000005","complex_substitution","SO:1000002");
		tree.newTerm("SO:1000008","point_mutation","SO:0001483");
		tree.newTerm("SO:1000009","transition","SO:0001483");
		tree.newTerm("SO:1000010","pyrimidine_transition","SO:1000009");
		tree.newTerm("SO:1000011","C_to_T_transition","SO:1000010");
		tree.newTerm("SO:1000012","C_to_T_transition_at_pCpG_site","SO:1000011");
		tree.newTerm("SO:1000013","T_to_C_transition","SO:1000010");
		tree.newTerm("SO:1000014","purine_transition","SO:1000009");
		tree.newTerm("SO:1000015","A_to_G_transition","SO:1000014");
		tree.newTerm("SO:1000016","G_to_A_transition","SO:1000014");
		tree.newTerm("SO:1000017","transversion","SO:0001483");
		tree.newTerm("SO:1000018","pyrimidine_to_purine_transversion","SO:1000017");
		tree.newTerm("SO:1000019","C_to_A_transversion","SO:1000018");
		tree.newTerm("SO:1000020","C_to_G_transversion","SO:1000018");
		tree.newTerm("SO:1000021","T_to_A_transversion","SO:1000018");
		tree.newTerm("SO:1000022","T_to_G_transversion","SO:1000018");
		tree.newTerm("SO:1000023","purine_to_pyrimidine_transversion","SO:1000017");
		tree.newTerm("SO:1000024","A_to_C_transversion","SO:1000023");
		tree.newTerm("SO:1000025","A_to_T_transversion","SO:1000023");
		tree.newTerm("SO:1000026","G_to_C_transversion","SO:1000023");
		tree.newTerm("SO:1000027","G_to_T_transversion","SO:1000023");
		tree.newTerm("SO:1000028","intrachromosomal_mutation","SO:1000183");
		tree.newTerm("SO:1000029","chromosomal_deletion","SO:1000028");
		tree.newTerm("SO:1000030","chromosomal_inversion","SO:1000028");
		tree.newTerm("SO:1000031","interchromosomal_mutation","SO:1000183");
		tree.newTerm("SO:1000032","indel","SO:0001059");
		tree.newTerm("SO:1000035","duplication","SO:0000667");
		tree.newTerm("SO:1000036","inversion","SO:0001059","SO:0001411");
		tree.newTerm("SO:1000037","chromosomal_duplication","SO:1000183");
		tree.newTerm("SO:1000038","intrachromosomal_duplication","SO:1000028","SO:1000037");
		tree.newTerm("SO:1000039","direct_tandem_duplication","SO:1000173");
		tree.newTerm("SO:1000040","inverted_tandem_duplication","SO:1000173");
		tree.newTerm("SO:1000041","intrachromosomal_transposition","SO:0000453","SO:1000038");
		tree.newTerm("SO:1000042","compound_chromosome","SO:1000183");
		tree.newTerm("SO:1000043","Robertsonian_fusion","SO:1000044");
		tree.newTerm("SO:1000044","chromosomal_translocation","SO:0000199","SO:1000031");
		tree.newTerm("SO:1000045","ring_chromosome","SO:1000028");
		tree.newTerm("SO:1000046","pericentric_inversion","SO:1000030");
		tree.newTerm("SO:1000047","paracentric_inversion","SO:1000030");
		tree.newTerm("SO:1000048","reciprocal_chromosomal_translocation","SO:1000044");
		tree.newTerm("SO:1000136","autosynaptic_chromosome","SO:1000183");
		tree.newTerm("SO:1000138","homo_compound_chromosome","SO:1000042");
		tree.newTerm("SO:1000140","hetero_compound_chromosome","SO:1000042");
		tree.newTerm("SO:1000141","chromosome_fission","SO:1000028");
		tree.newTerm("SO:1000142","dexstrosynaptic_chromosome","SO:1000136");
		tree.newTerm("SO:1000143","laevosynaptic_chromosome","SO:1000136");
		tree.newTerm("SO:1000144","free_duplication","SO:1000037");
		tree.newTerm("SO:1000145","free_ring_duplication","SO:1000045","SO:1000144");
		tree.newTerm("SO:1000147","deficient_translocation","SO:1000029","SO:1000044");
		tree.newTerm("SO:1000148","inversion_cum_translocation","SO:1000030","SO:1000044");
		tree.newTerm("SO:1000149","bipartite_duplication","SO:1000031","SO:1000038");
		tree.newTerm("SO:1000150","cyclic_translocation","SO:0002060");
		tree.newTerm("SO:1000151","bipartite_inversion","SO:1000030");
		tree.newTerm("SO:1000152","uninverted_insertional_duplication","SO:1000154");
		tree.newTerm("SO:1000153","inverted_insertional_duplication","SO:1000154");
		tree.newTerm("SO:1000154","insertional_duplication","SO:1000037");
		tree.newTerm("SO:1000155","interchromosomal_transposition","SO:0000453","SO:1000031");
		tree.newTerm("SO:1000156","inverted_interchromosomal_transposition","SO:1000155");
		tree.newTerm("SO:1000157","uninverted_interchromosomal_transposition","SO:1000155");
		tree.newTerm("SO:1000158","inverted_intrachromosomal_transposition","SO:1000148");
		tree.newTerm("SO:1000159","uninverted_intrachromosomal_transposition","SO:1000041");
		tree.newTerm("SO:1000160","unoriented_insertional_duplication","SO:1000154");
		tree.newTerm("SO:1000161","unoriented_interchromosomal_transposition","SO:1000155");
		tree.newTerm("SO:1000162","unoriented_intrachromosomal_transposition","SO:1000041");
		tree.newTerm("SO:1000170","uncharacterized_chromosomal_mutation","SO:1000183");
		tree.newTerm("SO:1000171","deficient_inversion","SO:1000029","SO:1000030");
		tree.newTerm("SO:1000173","tandem_duplication","SO:1000035");
		tree.newTerm("SO:1000175","partially_characterized_chromosomal_mutation","SO:1000170");
		tree.newTerm("SO:1000182","chromosome_number_variation","SO:0000240");
		tree.newTerm("SO:1000183","chromosome_structure_variation","SO:0000240");
		tree.newTerm("SO:1001187","alternatively_spliced_transcript","SO:0000673");
		tree.newTerm("SO:1001188","encodes_1_polypeptide","SO:0000463");
		tree.newTerm("SO:1001189","encodes_greater_than_1_polypeptide","SO:0000463");
		tree.newTerm("SO:1001190","encodes_different_polypeptides_different_stop","SO:1001195");
		tree.newTerm("SO:1001191","encodes_overlapping_peptides_different_start","SO:1001195");
		tree.newTerm("SO:1001192","encodes_disjoint_polypeptides","SO:1001189");
		tree.newTerm("SO:1001193","encodes_overlapping_polypeptides_different_start_and_stop","SO:1001195");
		tree.newTerm("SO:1001195","encodes_overlapping_peptides","SO:1001189");
		tree.newTerm("SO:1001196","cryptogene","SO:0000654","SO:0001431");
		tree.newTerm("SO:1001197","dicistronic_primary_transcript","SO:0000079","SO:0000631");
		tree.newTerm("SO:1001217","member_of_regulon","SO:0000081");
		tree.newTerm("SO:1001246","CDS_independently_known","SO:0000316");
		tree.newTerm("SO:1001247","orphan_CDS","SO:1001254");
		tree.newTerm("SO:1001249","CDS_supported_by_domain_match_data","SO:1001251");
		tree.newTerm("SO:1001251","CDS_supported_by_sequence_similarity_data","SO:1001254");
		tree.newTerm("SO:1001254","CDS_predicted","SO:0000316");
		tree.newTerm("SO:1001259","CDS_supported_by_EST_or_cDNA_data","SO:1001251");
		tree.newTerm("SO:1001260","internal_Shine_Dalgarno_sequence","SO:0000243","SO:1001268");
		tree.newTerm("SO:1001261","recoded_mRNA","SO:0000234");
		tree.newTerm("SO:1001262","minus_1_translationally_frameshifted","SO:0000887");
		tree.newTerm("SO:1001263","plus_1_translationally_frameshifted","SO:0000887");
		tree.newTerm("SO:1001264","mRNA_recoded_by_translational_bypass","SO:1001261");
		tree.newTerm("SO:1001265","mRNA_recoded_by_codon_redefinition","SO:1001261");
		tree.newTerm("SO:1001268","recoding_stimulatory_region","SO:0000836");
		tree.newTerm("SO:1001269","four_bp_start_codon","SO:0000680");
		tree.newTerm("SO:1001271","archaeal_intron","SO:0001216");
		tree.newTerm("SO:1001272","tRNA_intron","SO:0001216");
		tree.newTerm("SO:1001273","CTG_start_codon","SO:0000680");
		tree.newTerm("SO:1001274","SECIS_element","SO:1001268");
		tree.newTerm("SO:1001275","retron","SO:0001411");
		tree.newTerm("SO:1001277","three_prime_recoding_site","SO:1001268");
		tree.newTerm("SO:1001279","three_prime_stem_loop_structure","SO:1001277");
		tree.newTerm("SO:1001280","five_prime_recoding_site","SO:1001268");
		tree.newTerm("SO:1001281","flanking_three_prime_quadruplet_recoding_signal","SO:1001277");
		tree.newTerm("SO:1001282","UAG_stop_codon_signal","SO:1001288");
		tree.newTerm("SO:1001283","UAA_stop_codon_signal","SO:1001288");
		tree.newTerm("SO:1001284","regulon","SO:0005855");
		tree.newTerm("SO:1001285","UGA_stop_codon_signal","SO:1001288");
		tree.newTerm("SO:1001286","three_prime_repeat_recoding_signal","SO:1001277");
		tree.newTerm("SO:1001287","distant_three_prime_recoding_signal","SO:1001277");
		tree.newTerm("SO:1001288","stop_codon_signal","SO:1001268");
		tree.newTerm("SO:2000061","databank_entry","SO:0000695");
		tree.newTerm("SO:3000000","gene_segment","SO:0000842");		
		
		for(final TermImpl t:tree.acn2term.values()) {
			if(t.label==null) throw new JvarkitException.ProgrammingError("term "+t.accession+" has no label");
			if(t!=tree.getTermByLabel(t.label)) throw new JvarkitException.ProgrammingError("???");
		}
		return tree;
	 }
		
	
	private static class OwlLoader
		{
		private final String OWL="http://www.w3.org/2002/07/owl#";
		private final String OBOINOWL="http://www.geneontology.org/formats/oboInOwl#";
		private final String RDFS="http://www.w3.org/2000/01/rdf-schema#";
		private final String RDF="http://www.w3.org/1999/02/22-rdf-syntax-ns#";
		private SequenceOntologyTree tree = new SequenceOntologyTree();
		private final Map<String,TermImpl> uri2terms=new HashMap<>();
		private final QName rdfAbout=new QName(RDF, "about","rdf");
		private final QName rdfResource=new QName(RDF, "resource","rdf");
		private final List<String> uri2sub=new ArrayList<>();
		private void owlClass(final XMLEventReader r,final StartElement root) throws XMLStreamException {
			String acn=null;
			String label=null;
			final Attribute att = root.getAttributeByName(rdfAbout);
			if(att==null) throw new XMLStreamException("no rdf:about",root.getLocation());
			final String uri = att.getValue();
			boolean deprecated=false;
			while(r.hasNext()) {
				final XMLEvent evt= r.nextEvent();
				if( evt.isStartElement() ) {
					final StartElement E = evt.asStartElement();
					final QName qname = E.getName();
					if(qname.getLocalPart().equals("label") && RDFS.equals(qname.getNamespaceURI())) {
						label = r.getElementText();
						}
					else if(qname.getLocalPart().equals("id") && OBOINOWL.equals(qname.getNamespaceURI())) {
						acn = r.getElementText();
						}
					else if(qname.getLocalPart().equals("deprecated") && OWL.equals(qname.getNamespaceURI())) {
						if(r.getElementText().equals("true")) deprecated=true;
						}
					else if(qname.getLocalPart().equals("subClassOf") && RDFS.equals(qname.getNamespaceURI()) && E.getAttributeByName(rdfResource)!=null) {
						final Attribute att2 = E.getAttributeByName(rdfResource);
						if(att2==null) throw new XMLStreamException("no rdf:resource",E.getLocation());
						this.uri2sub.add(uri);
						this.uri2sub.add(att2.getValue());
						}
					}
				else if(evt.isEndElement()) {
					final EndElement E = evt.asEndElement();
					final QName qname = E.getName();
					if(qname.getLocalPart().equals("Class") && OWL.equals(qname.getNamespaceURI())) {
						break;
						}
					}
				}
			if(deprecated) return;
			if(acn==null)  throw new XMLStreamException("no acn",root.getLocation());
			if(label==null)  throw new XMLStreamException("no label",root.getLocation());
			if(this.uri2terms.containsKey(uri)) throw new XMLStreamException("duplicate uri ?" + uri,root.getLocation());
			final TermImpl term = tree.newTerm(acn,label);
			this.uri2terms.put(uri,term);
			}
		
		SequenceOntologyTree parse(final String uri) throws IOException {
			final XMLInputFactory xif = XMLInputFactory.newFactory();
			InputStream in = null;
			XMLEventReader r = null;
			try
				{
				in = IOUtils.openURIForReading(uri);
				r = xif.createXMLEventReader(in);
				return parse(r);
				}
			catch (final Exception e)
				{
				throw new RuntimeIOException(e);
				}
			finally  {
				CloserUtil.close(r);
				CloserUtil.close(in);
				}
			}
		
		SequenceOntologyTree parse(final InputStream in) throws RuntimeIOException {
			final XMLInputFactory xif = XMLInputFactory.newFactory();
			XMLEventReader r = null;
			try
				{
				r = xif.createXMLEventReader(in);
				return parse(r);
				}
			catch (final Exception e)
				{
				throw new RuntimeIOException(e);
				}
			finally  {
				CloserUtil.close(r);
				}
			}
		
		SequenceOntologyTree parse(final XMLEventReader r) throws XMLStreamException {
				this.tree = new SequenceOntologyTree();
				while(r.hasNext()) {
					final XMLEvent evt= r.nextEvent();
					if( evt.isStartElement() ) {
						final StartElement E = evt.asStartElement();
						final QName qname = E.getName();
						if(qname.getLocalPart().equals("Class") && OWL.equals(qname.getNamespaceURI())) {
							this.owlClass(r, E);
							}
						}
					}
				for(int i=0;i+1<this.uri2sub.size();i+=2) {
					final TermImpl child = this.uri2terms.get(this.uri2sub.get(i));
					if(child==null) throw new  XMLStreamException("cannot get child term "+ this.uri2sub.get(i));
					final TermImpl parent = this.uri2terms.get(this.uri2sub.get(i+1));
					if(parent==null) throw new  XMLStreamException("cannot get parent term "+ this.uri2sub.get(i+1));
					parent.children.add(child);
					child.parents.add(parent);
					}
				final SequenceOntologyTree t2 = this.tree;
				this.tree = null;
				this.uri2terms.clear();
				return t2;
				}
	
		}
	
	/** damaging comparator, try to order the term, the most damaging first */
	public static class DamagingComparator
		implements Comparator<Term>
		{
		private final Map<Term,Integer> term2score = new HashMap<>();
		
		private void initMap(final Term t,int score)
			{
			if(term2score.containsKey(Objects.requireNonNull(t, "term is  null"))) return;
			term2score.put(t, score);
			for(final Term c:t.getChildren()) {
				initMap(c,score+1);
				}
			}
		
		public DamagingComparator(final SequenceOntologyTree tree) {
			int startscore=1000;
			initMap(tree.getTermByLabel("missense_variant"),startscore);
			initMap(tree.getTermByLabel("nonsynonymous_variant"),startscore-=100);
			initMap(tree.getTermByLabel("protein_altering_variant"),startscore-=100);
			initMap(tree.getTermByLabel("coding_transcript_variant"),startscore-=100);
			initMap(tree.getTermByLabel("transcript_variant"),startscore-=100);
			
			initMap(tree.getTermByLabel("gene_variant"),startscore-=100);
			initMap(tree.getTermByLabel("regulatory_region_variant"),startscore);
			initMap(tree.getTermByLabel("intergenic_variant"),startscore);
			initMap(tree.getTermByLabel("intergenic_region"),startscore);
			//
			initMap(tree.getTermByLabel("sequence_variant"),100);
			initMap(tree.getTermByLabel("structural_variant"),300);
			}
		
		public  DamagingComparator() {
			this(SequenceOntologyTree.getInstance());
			}
		protected int computeDamagingScore(final Term t) {
			
			LOG.warn("Cannot compute damaging score of "+t);
			return 0;
			}
		/** the higher, the most precise/damaging */
		protected int getDamagingScore(final Term t)
			{
			Integer score = this.term2score.get(t);
			if(score==null) {
				score = computeDamagingScore(t);
				this.term2score.put(t,score);
				}
			return score.intValue();
			}
		@Override
		public int compare(final Term o1,final Term o2) {
			return Integer.compare(getDamagingScore(o2),getDamagingScore(o1));
			}
		}
	
	
	}
