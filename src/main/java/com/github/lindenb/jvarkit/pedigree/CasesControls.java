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
package com.github.lindenb.jvarkit.pedigree;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.CharSplitter;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.log.Logger;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.vcf.VCFHeader;


/**
 * A class parsing cases an controls
 */
public class CasesControls {
	private static final Logger LOG = Logger.of(CasesControls.class);

	@Parameter(names={"--cases"},description="File or comma-separated list of control samples")
	private String sourceCases = null;
	@Parameter(names={"--controls"},description="File or comma-separated list of control samples")
	private String sourceControls = null;

	private Set<String> _cases = Collections.emptySet();
	private Set<String> _controls = Collections.emptySet();
	
	
	public CasesControls() {
		}
	
	/** copy constructor */
	public CasesControls(final CasesControls cp) {
		this.sourceCases = cp.sourceCases;
		this.sourceControls = cp.sourceControls;
		this._cases=new HashSet<>(cp._cases);
		this._controls=new HashSet<>(cp._controls);
		}
	
	
	
	private Set<String> load(final String pathOrList) throws IOException {
		if(StringUtils.isBlank(pathOrList) ) return Collections.emptySet();
		final Path path = Paths.get(pathOrList);
		final Set<String> set;
		final String nature;
		if(Files.exists(path) && !Files.isDirectory(path) && Files.isReadable(path)) {
			nature = " (argument as file)";
			try(BufferedReader br = IOUtils.openPathForBufferedReading(path)) {
				set= br.lines().
						filter(S->!StringUtils.isBlank(S) || S.startsWith("#")).
						collect(Collectors.toSet())
						;
				}
			}
		else
			{
			nature = " (argument as comma-sparated)";
			set = Arrays.stream(CharSplitter.COMMA.split(pathOrList)).
					map(T->T.trim()).
					filter(S->!StringUtils.isBlank(S)).
					collect(Collectors.toSet())
					;
			}
		
		
		if(set.isEmpty()) LOG.warn("No valid sample was found in "+path+". "+nature);
		return Collections.unmodifiableSet(set);
		}
	public CasesControls load() {
		try {
			this._cases = load(this.sourceCases);
			this._controls = load(this.sourceControls);
			}
		catch(final IOException err) {
			throw new RuntimeIOException(err);
			}
		return this;
		}
	
	/** retain cases / controls that are in the VCF header */
	public CasesControls retain(final VCFHeader header) {
		if(header==null) throw new NullPointerException("vcf header is null");
		if(!header.hasGenotypingData()) {
			LOG.warn("the vcf header contains no genotype/sample.");
			}
		return retain(header.getGenotypeSamples());
		}
	
	/** retain cases / controls that are in collection */
	public CasesControls retain(final Collection<String> other) {
		if(other==null) throw new NullPointerException("other header is null");
		final Set<String> vcfsamples = new HashSet<>(other);
		if(vcfsamples.isEmpty()) {
			LOG.warn("CasesControls::retain the collection contains no genotype/sample.");
			}
		
		this._cases = Collections.unmodifiableSet(this._cases.stream().filter(S->vcfsamples.contains(S)).collect(Collectors.toSet()));
		this._controls = Collections.unmodifiableSet(this._controls.stream().filter(S->vcfsamples.contains(S)).collect(Collectors.toSet()));
		return this;
		}
	
	/** get Cases AND controls */
	public Set<String> getAll() {
		final Set<String> set = new TreeSet<>();
		set.addAll(getCases());
		set.addAll(getControls());
		return set;
		}

	
	public Set<String> getCases() {
		return _cases;
		}
	
	public int getCasesCount() {
		return getCases().size();
		}

	
	public Set<String> getControls() {
		return _controls;
		}
	
	public int getControlsCount() {
		return getControls().size();
		}

	
	/**
	 * return getCases() for side==0 && getControls() for side==1
	 */
	public Set<String> get(int side) {
		switch(side) {
		case 0: return this.getCases();
		case 1: return this.getControls();
		default: throw new ArrayIndexOutOfBoundsException(side);
		}
	}
	
	public boolean haveCases() {
		return !this._cases.isEmpty();
		}
	
	public boolean haveControls() {
		return !this._controls.isEmpty();
		}
	
	public boolean haveCasesAndControls() {
		return haveCases() && haveControls();
		}

	
	public CasesControls checkHaveCases() {
		if(!haveCases()) throw new IllegalArgumentException("No case(s) was found");
		return this;
		}
	public CasesControls checkHaveControls() {
		if(!haveControls()) throw new IllegalArgumentException("No control(s) was found");
		return this;
		}
	
	public CasesControls checkHaveCasesControls() {
		return checkHaveCases().checkHaveControls();
		}
	
	/** return true if cases or controls contains(gt.getSampleName()) */
	public boolean contains(final Genotype gt) {
		return contains(gt.getSampleName());
	}
	
	/** return true if cases or controls contains sn */
	public boolean contains(final String sample) {
		return isCase(sample) || isControl(sample);
	}
	
	public boolean isCase(final Genotype gt) {
		return isCase(gt.getSampleName());
		}
	
	public boolean isControl(final Genotype gt) {
		return isControl(gt.getSampleName());
		}
	
	public boolean isCase(final String sn) {
		return this._cases.contains(sn);
		}
	
	public boolean isControl(final String sn) {
		return this._controls.contains(sn);
		}
	
	public CasesControls checkNoCommon() {
		final Set<String> common = _cases.stream().filter(S->_controls.contains(S)).collect(Collectors.toSet());
		if(_cases.stream().anyMatch(S->_controls.contains(S))) {
			throw new IllegalArgumentException("Common samples between cases and controls:"+String.join(" ", common));
			}
		return this;
		}
	
	/** return true if no case and no control */
	public boolean isEmpty() {
		return getCases().isEmpty() && getControls().isEmpty();
	}

	@Override
	public CasesControls clone() {
		return new CasesControls(this);
		}
	
	
	@Override
	public boolean equals(final Object obj) {
		if(obj==this) return true;
		if(obj==null || !(obj instanceof CasesControls)) return false;
		final CasesControls o = CasesControls.class.cast(obj);
		return o._cases.equals(this._cases) && o._controls.equals(this._controls);
		}
	
	@Override
	public int hashCode() {
		return _cases.hashCode()*31+ this._controls.hashCode();
		}
	
	@Override
	public String toString() {
		return "cases N="+ this._cases.size()+" : "+ String.join(" ", this._cases) +"\ncontrols N="+this._controls.size()+" : "+String.join(" ", this._controls);
		}
	}
