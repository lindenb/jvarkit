package com.github.lindenb.jvarkit.tools.cloneit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.util.log.Logger;

public class SubCloneIt extends AbstractCloneIt {
	private static final Logger LOG = Logger.build(SubCloneIt.class).make();
	
	@ParametersDelegate
	private InsertPlasmidLoader insertLoader=new InsertPlasmidLoader();
	@ParametersDelegate
	private VectorPlasmidLoader vectorLoader=new VectorPlasmidLoader();
	
	private float enzymeWeight=5f;
	private InsertPlasmid insert=null;
	private VectorPlasmid vector=null;
	private boolean acceptPolymerase=false;
	private boolean useCIAP=false;
	
private class HemiStrategy
	{
	SiteWithTreatment left;
	SiteWithTreatment right;
	}

private class Strategy
	{
	HemiStrategy left;
	HemiStrategy riht;
	}


@Override
public int doWork(final List<String> args) {
	if(!args.isEmpty()) {
		LOG.error("too many arguments");
		return -1;
	}
	try 
		{
		final List<Enzyme> enzymes = getRebase().stream().
				filter(E->E.isPalindromic()).
				filter(E->E.getWeight()>this.enzymeWeight).
				map(E->new EnzymeImpl(E)).
				collect(Collectors.toList());
		if(enzymes.isEmpty()) {
			LOG.error("enzyme list is empty");
			return -1;
		}
		this.insert  = this.insertLoader.create();
		this.vector  = this.vectorLoader.create();
		this.insert.digest(enzymes);
		this.vector.digest(enzymes);

		
		final Set<Enzyme> badEnzymesV = this.vector.getSites().stream().
				filter(S->S.getPosition()< this.vector.getPolylinker().getStart() || S.getPosition()> this.vector.getPolylinker().getEnd()).
				map(S->S.getEnzyme()).
				collect(Collectors.toSet())
				;
		
		final Set<Enzyme> badEnzymesI = this.insert.getSites().stream().
				filter(S->
						S.getPosition()< this.insert.getPolylinker5().getStart() ||
						S.getPosition()> this.insert.getPolylinker3().getEnd() ||
						(S.getPosition()> this.insert.getPolylinker5().getEnd()  && S.getPosition()< this.insert.getPolylinker3().getStart())
						).
				map(S->S.getEnzyme()).
				collect(Collectors.toSet())
				;
		
		
		
		final List<Site> listV= this.vector.getSites().stream().
				filter(S->S.getPosition()>=this.vector.getPolylinker().getStart()).
				filter(S->S.getPosition()<=this.vector.getPolylinker().getEnd()).
				filter(S->!badEnzymesV.contains(S.getEnzyme())).
				collect(Collectors.toList());
		
		if(listV.isEmpty()) {
			LOG.error("Can find any suitable enzyme in vector "+vector.getName());
			return -1;
		}
		
		final List<Site> listI5= this.insert.getSites().stream().
				filter(S->S.getPosition()>=this.insert.getPolylinker5().getStart()).
				filter(S->S.getPosition()<=this.insert.getPolylinker5().getEnd()).
				filter(S->!badEnzymesI.contains(S.getEnzyme())).
				collect(Collectors.toList());
		
		if(listI5.isEmpty()) {
			LOG.error("Can find any suitable enzyme in insert 5' "+insert.getName());
			return -1;
			}
		
		final List<Site> listI3= this.insert.getSites().stream().
				filter(S->S.getPosition()>=this.insert.getPolylinker3().getStart()).
				filter(S->S.getPosition()<=this.insert.getPolylinker3().getEnd()).
				filter(S->!badEnzymesI.contains(S.getEnzyme())).
				collect(Collectors.toList());
		
		if(listI3.isEmpty()) {
			LOG.error("Can find any suitable enzyme in insert 3' "+insert.getName());
			return -1;
			}
		
		listV.stream().forEach(siteV5->{
		listV.stream().
			filter(S->S.equals(siteV5) || S.getPosition()>siteV5.getMaxPosition()).
			forEach(siteV3->{
				
				
				
			});
			
			
		});
		
		for(final Site siteV5: listV) {
			
			for(final Site siteV3: listV) {
				
				for(final Site siteI5: listI5) {
					
					for(final Site siteI3: listI3) {
						
						Strategy strategy = new Strategy();
						
						
						}
					
					}

				
				}
			
			}
		
		
		return 0;
		}
	catch(final Exception err) {
		return -1;
		}
	
	}

public static void main(String[] args) {
	new SubCloneIt().instanceMainWithExit(args);
	}
}
