/*
The MIT License (MIT)
Copyright (c) 2024 Pierre Lindenbaum

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
package com.github.lindenb.jvarkit.tools.gatk;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import com.github.lindenb.jvarkit.gatk.Gatk4Proxy;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.util.FileExtensions;

/**
 * Base class for tools using the GATK wrapper
 * @author lindenb
 *
 */
public class AbstractGatkTool extends Launcher {
	private static final Logger LOG = Logger.build(AbstractGatkTool.class).make();
	private final Gatk4Proxy gatkEngine = Gatk4Proxy.getInstance().orElse(null);
	
	protected Gatk4Proxy getGatkEngine() {
		return gatkEngine;
		}
	
	protected void execute(final List<String> argv) throws Exception {
		getGatkEngine().execute(argv);
		}
	/** retrieve path for index of given hts file (vcf.gz, bam, cram...) */
	protected Path indexFor(final Path base) {
		final String fname = base.getFileName().toString();
		if(fname.endsWith(FileExtensions.COMPRESSED_VCF)) {
			return base.getParent().resolve(fname+FileExtensions.TABIX_INDEX);
			}
		else if(fname.endsWith(FileExtensions.VCF)) {
			return base.getParent().resolve(fname+FileExtensions.TRIBBLE_INDEX);
			}
		else if(fname.endsWith(FileExtensions.BAM)) {
			return base.getParent().resolve(fname+FileExtensions.BAI_INDEX);
			}
		else if(fname.endsWith(FileExtensions.CRAM)) {
			return base.getParent().resolve(fname+FileExtensions.CRAM_INDEX);
			}
		throw new IllegalArgumentException("Cannot find index for "+base);
		}
	
	protected void FilesDelete(final Path path) throws IOException {
		getLogger().info("delete "+path);
		Files.delete(path);
		}
	protected void FilesDeleteIfExists(final Path path) throws IOException {
		getLogger().info("delete "+path);
		Files.deleteIfExists(path);
		}
	
	protected Logger getLogger() {
		return LOG;
		}
	}
