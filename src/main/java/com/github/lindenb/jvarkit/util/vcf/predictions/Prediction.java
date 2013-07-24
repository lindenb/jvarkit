package com.github.lindenb.jvarkit.util.vcf.predictions;

import java.util.Set;

import com.github.lindenb.jvarkit.util.so.SequenceOntologyTree;

public interface Prediction {
public String getGeneName();
public String getEnsemblGene();
public String getEnsemblTranscript();
public String getEnsemblProtein();
public String getReferenceAminoAcid();
public Integer getAminoAcidPosition();
public String getAltAminoAcid();
public Set<SequenceOntologyTree.Term> getSOTerms();
}
