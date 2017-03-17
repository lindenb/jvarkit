/*
The MIT License (MIT)

Copyright (c) 2014 Pierre Lindenbaum

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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import com.github.lindenb.jvarkit.io.IOUtils;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;



/**Sequence ontology tree */
public class SequenceOntologyTree
	{

	private static SequenceOntologyTree INSTANCE=null;
	private final Map<String,TermImpl> uri2term=new HashMap<String,TermImpl>();
	private final Map<String,TermImpl> label2term=new HashMap<String,TermImpl>();
	
	
	public interface Term
		{
		public String getUri();
		public String getAcn();
		public String getLabel();
		public Set<Term> getParents();
		/** get direct children of this node */
		public Set<Term> getChildren();
		/** get ALL (recursive) children of this node */
		public Set<Term> getAllDescendants();

		}
	private class TermImpl implements Term
		{
		final String accession;
		final String name;
		final Set<String> parents=new HashSet<String>();
		final Set<String> children=new HashSet<String>();
		
		TermImpl(final String accession,final String name) {
			this.accession = accession;
			this.name = name;
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
			return name;
			}
		
		private Set<Term> convert(final Set<String> S1)
			{
			final Set<Term> S2=new HashSet<Term>(S1.size());
			for(final String s:S1)
				{
				final Term t=uri2term.get(s);
				if(t==null) continue;
				S2.add(t);
				}
			return S2;
			}
		
		/** returns only the direct children of this node */
		@Override
		public Set<Term> getChildren()
			{
			return convert(children);
			}
		
		/** recursive operation on getChildren, including self */
		@Override
		public Set<Term> getAllDescendants()
			{
			final Set<Term> set=new HashSet<Term>();
			_getAllDescendants(this,set);
			return set;
			}
		
		
		@Override
		public Set<Term> getParents()
			{
			return convert(parents);
			}
		
		@Override
		public int hashCode()
			{
			return accession.hashCode();
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
	
	private static void _getAllDescendants(Term term,Set<Term> set)
		{
		set.add(term);
		for(Term c:term.getChildren())
			{
			_getAllDescendants(c,set);
			}
		}

	private TermImpl createTerm(final String acn,final String label)
		{
		final TermImpl t=new TermImpl(acn,label);
		this.uri2term.put(t.accession, t);
		this.label2term.put(t.getLabel(), t);
		return t;
		}
	
	private void addTerm(final String acn,final String label,final String parent,final String children)
		{
		final TermImpl t= createTerm(acn,label);
		final Pattern comma=Pattern.compile("[,]");
		for(final String s:comma.split(parent))
			{
			if(s.isEmpty()) continue;
			t.parents.add(s);
			}
		for(final String s:comma.split(children))
			{
			if(s.isEmpty()) continue;
			t.children.add(s);
			}
		}
	
	public Term getTermByAcn(final String s)
		{
		return this.uri2term.get(s);
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
		return this.uri2term.values();
		}
	
	public String getLastUpdated()
		{
		return "2015-01-23";
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
					INSTANCE= createDefault();
					}
				}
			}
		return INSTANCE;
		}
	public static SequenceOntologyTree createFromInputStream(final InputStream in) throws IOException
		{
		return new OwlLoader().parse(in);
		}
	public static SequenceOntologyTree createDefault()
		{
		final SequenceOntologyTree tree =new SequenceOntologyTree();
						//INSTANCE.addTerm("Accession","Name","DB_Xrefs","Parents","Children");
		tree.addTerm("SO:0001628","intergenic_variant","SO:0001878","SO:0001632,SO:0001986,SO:0002017,SO:0001631,SO:0001987");
		tree.addTerm("SO:0001782","TF_binding_site_variant","SO:0001566","");
		tree.addTerm("SO:0001565","gene_fusion","SO:0001564","");
		tree.addTerm("SO:0001826","disruptive_inframe_deletion","SO:0001822","");
		tree.addTerm("SO:0001537","structural_variant","SO:0001060","SO:0001880,SO:0001882,SO:0001879,SO:0001878,SO:0001881,SO:0001563");
		tree.addTerm("SO:0001591","frame_restoring_variant","SO:0001589","");
		tree.addTerm("SO:0001610","elongated_polypeptide_C_terminal","SO:0001609","SO:0001612,SO:0001613");
		tree.addTerm("SO:0001912","copy_number_decrease","SO:0001563","");
		tree.addTerm("SO:0002012","start_lost","SO:0001582,SO:0001992","");
		tree.addTerm("SO:0001060","sequence_variant","","SO:0001536,SO:0001537");
		tree.addTerm("SO:0001575","splice_donor_variant","SO:0001629","");
		tree.addTerm("SO:0001541","decreased_transcript_level_variant","SO:0001540","");
		tree.addTerm("SO:0001583","missense_variant","SO:0001992","SO:0001586,SO:0001585");
		tree.addTerm("SO:0001560","inactive_ligand_binding_site","SO:0001559","SO:0001618");
		tree.addTerm("SO:0001552","decreased_transcription_rate_variant","SO:0001550","");
		tree.addTerm("SO:0001787","splice_donor_5th_base_variant","SO:0001629","");
		tree.addTerm("SO:0001555","decreased_translational_product_level","SO:0001553","");
		tree.addTerm("SO:0001818","protein_altering_variant","SO:0001580","SO:0001650,SO:0001589");
		tree.addTerm("SO:0001791","exon_variant","SO:0001576","SO:0001580,SO:0001792,SO:0001622");
		tree.addTerm("SO:0001550","rate_of_transcription_variant","SO:0001549","SO:0001551,SO:0001552");
		tree.addTerm("SO:0001894","regulatory_region_ablation","SO:0001879","SO:0001895");
		tree.addTerm("SO:0001992","nonsynonymous_variant","SO:0001650","SO:0001583,SO:0001587,SO:0002012,SO:0001578");
		tree.addTerm("SO:0001822","inframe_deletion","SO:0001820,SO:0001906","SO:0001826,SO:0001825");
		tree.addTerm("SO:0001803","decreased_polyadenylation_variant","SO:0001545","");
		tree.addTerm("SO:0001824","disruptive_inframe_insertion","SO:0001821","");
		tree.addTerm("SO:0001605","amino_acid_insertion","SO:0001603","");
		tree.addTerm("SO:0001626","incomplete_terminal_codon_variant","SO:0001650,SO:0001590","");
		tree.addTerm("SO:0001632","downstream_gene_variant","SO:0001628","SO:0001633");
		tree.addTerm("SO:0001574","splice_acceptor_variant","SO:0001629","");
		tree.addTerm("SO:0001585","conservative_missense_variant","SO:0001583","");
		tree.addTerm("SO:0001895","TFBS_ablation","SO:0001894","");
		tree.addTerm("SO:0001612","elongated_in_frame_polypeptide_C_terminal","SO:0001610","");
		tree.addTerm("SO:0001616","polypeptide_fusion","SO:0001603","");
		tree.addTerm("SO:0001882","feature_fusion","SO:0001537","SO:0001886,SO:0001887,SO:0001890");
		tree.addTerm("SO:0001547","decreased_transcript_stability_variant","SO:0001546","");
		tree.addTerm("SO:0001593","minus_2_frameshift_variant","SO:0001589","");
		tree.addTerm("SO:0001578","stop_lost","SO:0001907,SO:0001590,SO:0001992","");
		tree.addTerm("SO:0001598","translational_product_structure_variant","SO:0001564","SO:0001603");
		tree.addTerm("SO:0001990","five_prime_UTR_premature_start_codon_location_variant","SO:0001983","");
		tree.addTerm("SO:0002010","pyrrolysine_loss","SO:0002008","");
		tree.addTerm("SO:0001613","elongated_out_of_frame_polypeptide_C_terminal","SO:0001610","");
		tree.addTerm("SO:0002014","5_prime_UTR_elongation","SO:0001623","");
		tree.addTerm("SO:0001878","feature_variant","SO:0001537","SO:0001628,SO:0001564,SO:0001907,SO:0001566,SO:0001017,SO:0001906");
		tree.addTerm("SO:0001562","polypeptide_post_translational_processing_variant","SO:0001554","");
		tree.addTerm("SO:0001987","downstream_transcript_variant","SO:0001628","");
		tree.addTerm("SO:0001908","internal_feature_elongation","SO:0001907","SO:0001909,SO:0001821");
		tree.addTerm("SO:0001970","non_coding_transcript_intron_variant","SO:0001619,SO:0001627","");
		tree.addTerm("SO:0001621","NMD_transcript_variant","SO:0001576","");
		tree.addTerm("SO:0001587","stop_gained","SO:0001992,SO:0001906","");
		tree.addTerm("SO:0001615","elongated_out_of_frame_polypeptide_N_terminal","SO:0001611","");
		tree.addTerm("SO:0001542","increased_transcript_level_variant","SO:0001540","");
		tree.addTerm("SO:0001636","2KB_upstream_variant","SO:0001635","");
		tree.addTerm("SO:0001968","coding_transcript_variant","SO:0001576","SO:0001969,SO:0001580,SO:0001622");
		tree.addTerm("SO:0001548","increased_transcript_stability_variant","SO:0001546","");
		tree.addTerm("SO:0001609","elongated_polypeptide","SO:0001603","SO:0001610,SO:0001611");
		tree.addTerm("SO:0001559","polypeptide_loss_of_function_variant","SO:0001554","SO:0001561,SO:0001560");
		tree.addTerm("SO:0001557","polypeptide_gain_of_function_variant","SO:0001554","");
		tree.addTerm("SO:0001546","transcript_stability_variant","SO:0001538","SO:0001547,SO:0001548");
		tree.addTerm("SO:0001989","5_prime_UTR_premature_start_codon_loss_variant","SO:0001983","");
		tree.addTerm("SO:0002008","rare_amino_acid_variant","SO:0001586","SO:0002009,SO:0002010");
		tree.addTerm("SO:0001619","non_coding_transcript_variant","SO:0001576","SO:0001620,SO:0001792,SO:0001970");
		tree.addTerm("SO:0001549","transcription_variant","SO:0001538","SO:0001550");
		tree.addTerm("SO:0001884","regulatory_region_translocation","SO:0001881","SO:0001885");
		tree.addTerm("SO:0001792","non_coding_transcript_exon_variant","SO:0001619,SO:0001791","");
		tree.addTerm("SO:0001544","editing_variant","SO:0001543","");
		tree.addTerm("SO:0001889","transcript_amplification","SO:0001880","");
		tree.addTerm("SO:0001880","feature_amplification","SO:0001537","SO:0001891,SO:0001889");
		tree.addTerm("SO:0001571","cryptic_splice_donor","SO:0001569","");
		tree.addTerm("SO:0001892","TFBS_amplification","SO:0001891","");
		tree.addTerm("SO:0001802","increased_polyadenylation_variant","SO:0001545","");
		tree.addTerm("SO:0001554","polypeptide_function_variant","SO:0001539","SO:0001557,SO:0001562,SO:0001558,SO:0001559");
		tree.addTerm("SO:0001909","frameshift_elongation","SO:0001589,SO:0001908","");
		tree.addTerm("SO:0001886","transcript_fusion","SO:0001882","");
		tree.addTerm("SO:0001595","plus_2_frameshift_variant","SO:0001589","");
		tree.addTerm("SO:0001539","translational_product_function_variant","SO:0001536","SO:0001602,SO:0001553,SO:0001599,SO:0001554");
		tree.addTerm("SO:0001879","feature_ablation","SO:0001537","SO:0001893,SO:0001894");
		tree.addTerm("SO:0001589","frameshift_variant","SO:0001818","SO:0001594,SO:0001592,SO:0001591,SO:0001909,SO:0001595,SO:0001593,SO:0001910");
		tree.addTerm("SO:0002019","start_retained_variant","SO:0001819,SO:0001582","");
		tree.addTerm("SO:0001906","feature_truncation","SO:0001878","SO:0001587,SO:0001910,SO:0001822");
		tree.addTerm("SO:0001567","stop_retained_variant","SO:0001819,SO:0001590","");
		tree.addTerm("SO:0001597","compensatory_transcript_secondary_structure_variant","SO:0001596","");
		tree.addTerm("SO:0001627","intron_variant","SO:0001576","SO:0001969,SO:0001629,SO:0002018,SO:0001970");
		tree.addTerm("SO:0001911","copy_number_increase","SO:0001563","");
		tree.addTerm("SO:0001604","amino_acid_deletion","SO:0001603","");
		tree.addTerm("SO:0001540","level_of_transcript_variant","SO:0001538","SO:0001542,SO:0001541");
		tree.addTerm("SO:0001551","increased_transcription_rate_variant","SO:0001550","");
		tree.addTerm("SO:0001592","minus_1_frameshift_variant","SO:0001589","");
		tree.addTerm("SO:0001573","intron_gain","SO:0001568","");
		tree.addTerm("SO:0002011","intragenic_variant","SO:0001576","");
		tree.addTerm("SO:0001611","elongated_polypeptide_N_terminal","SO:0001609","SO:0001614,SO:0001615");
		tree.addTerm("SO:0002018","conserved_intron_variant","SO:0001627","");
		tree.addTerm("SO:0001630","splice_region_variant","SO:0001568","");
		tree.addTerm("SO:0001553","translational_product_level_variant","SO:0001539","SO:0001555,SO:0001556");
		tree.addTerm("SO:0001576","transcript_variant","SO:0001564","SO:0001621,SO:0001627,SO:0002011,SO:0001619,SO:0001791,SO:0001968,SO:0001577,SO:0001596,SO:0001568");
		tree.addTerm("SO:0001819","synonymous_variant","SO:0001580","SO:0002019,SO:0001567");
		tree.addTerm("SO:0001538","transcript_function_variant","SO:0001536","SO:0001543,SO:0001540,SO:0001549,SO:0001546");
		tree.addTerm("SO:0001590","terminator_codon_variant","SO:0001580","SO:0001626,SO:0001578,SO:0001567");
		tree.addTerm("SO:0001601","conformational_change_variant","SO:0001599","");
		tree.addTerm("SO:0001545","polyadenylation_variant","SO:0001543","SO:0001803,SO:0001802");
		tree.addTerm("SO:0001603","polypeptide_sequence_variant","SO:0001598","SO:0001604,SO:0001617,SO:0001616,SO:0001605,SO:0001609,SO:0001606");
		tree.addTerm("SO:0001570","cryptic_splice_acceptor","SO:0001569","");
		tree.addTerm("SO:0002013","5_prime_UTR_truncation","SO:0001623","");
		tree.addTerm("SO:0001986","upstream_transcript_variant","SO:0001628","");
		tree.addTerm("SO:0001607","conservative_amino_acid_substitution","SO:0001606","");
		tree.addTerm("SO:0001582","initiator_codon_variant","SO:0001580","SO:0002019,SO:0002012");
		tree.addTerm("SO:0001983","5_prime_UTR_premature_start_codon_variant","SO:0001623","SO:0001988,SO:0001989,SO:0001990");
		tree.addTerm("SO:0001988","5_prime_UTR_premature_start_codon_gain_variant","SO:0001983","");
		tree.addTerm("SO:0002016","3_prime_UTR_elongation","SO:0001624","");
		tree.addTerm("SO:0001825","conservative_inframe_deletion","SO:0001822","");
		tree.addTerm("SO:0001618","inactive_catalytic_site","SO:0001560","");
		tree.addTerm("SO:0001890","transcript_regulatory_region_fusion","SO:0001882","");
		tree.addTerm("SO:0001617","polypeptide_truncation","SO:0001603","");
		tree.addTerm("SO:0001629","splice_site_variant","SO:0001627,SO:0001568","SO:0001574,SO:0001575,SO:0001787");
		tree.addTerm("SO:0001569","cryptic_splice_site_variant","SO:0001568","SO:0001570,SO:0001571");
		tree.addTerm("SO:0001820","inframe_indel","SO:0001650","SO:0001821,SO:0001822");
		tree.addTerm("SO:0001017","silent_mutation","SO:0001878","");
		tree.addTerm("SO:0001561","polypeptide_partial_loss_of_function","SO:0001559","");
		tree.addTerm("SO:0001566","regulatory_region_variant","SO:0001878","SO:0001782");
		tree.addTerm("SO:0001580","coding_sequence_variant","SO:0001791,SO:0001968","SO:0001819,SO:0001818,SO:0001590,SO:0001582");
		tree.addTerm("SO:0001634","500B_downstream_variant","SO:0001633","");
		tree.addTerm("SO:0001786","loss_of_heterozygosity","SO:0001536","");
		tree.addTerm("SO:0001907","feature_elongation","SO:0001878","SO:0001578,SO:0001908");
		tree.addTerm("SO:0001623","5_prime_UTR_variant","SO:0001622","SO:0001983,SO:0002014,SO:0002013");
		tree.addTerm("SO:0001577","complex_transcript_variant","SO:0001576","");
		tree.addTerm("SO:0001568","splicing_variant","SO:0001576","SO:0001573,SO:0001995,SO:0001572,SO:0001629,SO:0001569,SO:0001630");
		tree.addTerm("SO:0001821","inframe_insertion","SO:0001820,SO:0001908","SO:0001823,SO:0001824");
		tree.addTerm("SO:0001564","gene_variant","SO:0001878","SO:0001576,SO:0001565,SO:0001598");
		tree.addTerm("SO:0001883","transcript_translocation","SO:0001881","");
		tree.addTerm("SO:0001558","polypeptide_localization_variant","SO:0001554","");
		tree.addTerm("SO:0001622","UTR_variant","SO:0001791,SO:0001968","SO:0001623,SO:0001624");
		tree.addTerm("SO:0001536","functional_variant","SO:0001060","SO:0001786,SO:0001539,SO:0001538");
		tree.addTerm("SO:0001881","feature_translocation","SO:0001537","SO:0001884,SO:0001883");
		tree.addTerm("SO:0001586","non_conservative_missense_variant","SO:0001583","SO:0002008");
		tree.addTerm("SO:0001543","transcript_processing_variant","SO:0001538","SO:0001545,SO:0001544");
		tree.addTerm("SO:0001823","conservative_inframe_insertion","SO:0001821","");
		tree.addTerm("SO:0001624","3_prime_UTR_variant","SO:0001622","SO:0002016,SO:0002015");
		tree.addTerm("SO:0001887","regulatory_region_fusion","SO:0001882","SO:0001888");
		tree.addTerm("SO:0001891","regulatory_region_amplification","SO:0001880","SO:0001892");
		tree.addTerm("SO:0001599","3D_polypeptide_structure_variant","SO:0001539","SO:0001601,SO:0001600");
		tree.addTerm("SO:0002009","selenocysteine_loss","SO:0002008","");
		tree.addTerm("SO:0001608","non_conservative_amino_acid_substitution","SO:0001606","");
		tree.addTerm("SO:0001888","TFBS_fusion","SO:0001887","");
		tree.addTerm("SO:0001620","mature_miRNA_variant","SO:0001619","");
		tree.addTerm("SO:0001606","amino_acid_substitution","SO:0001603","SO:0001608,SO:0001607");
		tree.addTerm("SO:0001633","5KB_downstream_variant","SO:0001632","SO:0001634");
		tree.addTerm("SO:0001594","plus_1_frameshift_variant","SO:0001589","");
		tree.addTerm("SO:0001910","frameshift_truncation","SO:0001589,SO:0001906","");
		tree.addTerm("SO:0001650","inframe_variant","SO:0001818","SO:0001626,SO:0001992,SO:0001820");
		tree.addTerm("SO:0002017","conserved_intergenic_variant","SO:0001628","");
		tree.addTerm("SO:0001635","5KB_upstream_variant","SO:0001631","SO:0001636");
		tree.addTerm("SO:0001572","exon_loss_variant","SO:0001568","");
		tree.addTerm("SO:0001995","extended_intronic_splice_region_variant","SO:0001568","");
		tree.addTerm("SO:0001631","upstream_gene_variant","SO:0001628","SO:0001635");
		tree.addTerm("SO:0001563","copy_number_change","SO:0001537","SO:0001911,SO:0001912");
		tree.addTerm("SO:0001600","complex_3D_structural_variant","SO:0001599","");
		tree.addTerm("SO:0002015","3_prime_UTR_truncation","SO:0001624","");
		tree.addTerm("SO:0001969","coding_transcript_intron_variant","SO:0001968,SO:0001627","");
		tree.addTerm("SO:0001893","transcript_ablation","SO:0001879","");
		tree.addTerm("SO:0001602","complex_change_of_translational_product_variant","SO:0001539","");
		tree.addTerm("SO:0001614","elongated_in_frame_polypeptide_N_terminal_elongation","SO:0001611","");
		tree.addTerm("SO:0001556","increased_translational_product_level","SO:0001553","");
		tree.addTerm("SO:0001596","transcript_secondary_structure_variant","SO:0001576","SO:0001597");
		tree.addTerm("SO:0001885","TFBS_translocation","SO:0001884","");
		//2017 March 17
		tree.addTerm("SO:0000605","intergenic_region", "", "");
		tree.addTerm("SO:0001792","non_coding_exon_variant", "SO:0001791", "");
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
			final TermImpl term = tree.createTerm(acn,label);;
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
					parent.children.add(child.accession);
					child.parents.add(parent.accession);
					}
				final SequenceOntologyTree t2 = this.tree;
				this.tree = null;
				this.uri2terms.clear();
				return t2;
				}
	
		}
	
	
	
	}
