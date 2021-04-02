/**
 * 
 */
package com.github.lindenb.jvarkit.jcommander;

import java.io.File;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParametersDelegate;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.io.NullOuputStream;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.util.JVarkitVersion;
import com.github.lindenb.jvarkit.util.jcommander.Launcher.WritingBamArgs;
import com.github.lindenb.jvarkit.util.jcommander.Launcher.WritingSamReaderType;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;

@SuppressWarnings("unused")
public abstract class AbstractBamSplitter<T> extends MultiBamLauncher {

	private static final Logger LOG = Logger.build(AbstractBamSplitter.class).make();

	@Parameter(names= {"-o","--output"},description="(prefix) output directory",required=true)
	private Path outputDir=null;
	@Parameter(names= {"-M","--manifest"},description="Manifest file describing the generated files. Optional")
	private Path manifestFile=null;
	@Parameter(names= {"--prefix"},description="Output file prefix")
	private String prefix="split";
	@Parameter(names= {"--force"},description="overwrite existing files")
	private boolean force_overwrite = false;
	@ParametersDelegate
	private WritingBamArgs writingBamArgs=new WritingBamArgs().setSamOutputFormat(WritingSamReaderType.BAM);
	
	private static class SAMCount {
		final SAMFileWriter sfw;
		long count  = 0L;
		Path path;
		SAMCount(final SAMFileWriter sfw) {
			this.sfw = sfw;
			}
		}
	
	@Override
	protected Logger getLogger() {
		return LOG;
	};
	
	@Override
	protected int beforeSam() {
		IOUtil.assertDirectoryIsWritable(this.outputDir);
		if(this.prefix.startsWith(".") || this.prefix.contains(File.separator)) {
			LOG.error("bad file prefix " + this.prefix);
			}
		return super.beforeSam();
		}
	
	protected abstract Set<T> createKeys(final SAMRecord rec);
	
	/** initializing things with incoming SAM header */
	protected int beforeIterator(final SAMFileHeader header) {
		return 0;
	}
	
	@Override
	protected int processInput(final SAMFileHeader header, final CloseableIterator<SAMRecord> iter) {
		if(beforeIterator(header)!=0) {
			getLogger().error("initialisation failed.");
			return -1;
		}
		
		this.writingBamArgs.setReferencePath(super.faidxPath);
		
		
		
		final Set<String> samples = header.getReadGroups().stream().map(RG->RG.getSample()).filter(S->!StringUtils.isBlank(S)).collect(Collectors.toSet());
		final String samplePrefix;
		if(samples.size()==1) {
			samplePrefix= "."+ samples.iterator().next();
		} else
		{
			samplePrefix= "";
		}
		
		final Map<T,SAMCount> tag2sw = new HashMap<>();
		try(PrintWriter mOut = this.manifestFile==null?new PrintWriter(new NullOuputStream()):IOUtils.openPathForPrintWriter(this.manifestFile))
			{
			long n_lost=0L;
			while(iter.hasNext()) {
				final SAMRecord rec = iter.next();
				final Set<T> keys = createKeys(rec);
				if(keys==null || keys.isEmpty()) {
					n_lost++;
					continue;
				}
				
				for(final T id : keys) {
					SAMCount sc = tag2sw.get(id);
					if(sc==null) {
						getLogger().info("Creating output for \""+id+"\" N="+(tag2sw.size()+1));
						final Path path = this.outputDir.resolve(this.prefix+ samplePrefix + String.format(".%06d",(tag2sw.size()+1))+this.writingBamArgs.getFileExtension());
						if(Files.exists(path) && !force_overwrite) {
							getLogger().error("file exists : "+ path+ ". Use --force to overwrite");
							return -1;
							}
						final SAMFileHeader h2 = header.clone();
						h2.addComment(getProgramName()+" : " + id.toString());
						JVarkitVersion.getInstance().addMetaData(this, h2);
						sc = new SAMCount(this.writingBamArgs.openSamWriter(path, h2, true));
						sc.path = path;
						tag2sw.put(id, sc);
						}
					sc.sfw.addAlignment(rec);
					sc.count++;
					}
				}
			if(n_lost>0) getLogger().warn(StringUtils.niceInt(n_lost)+" read(s) where lost because no group was assigned.");
			
			mOut.println("#KEY\tPATH\tCOUNT");
			for(final T id : tag2sw.keySet()) {
				final SAMCount sc = tag2sw.get(id);
				mOut.print(id.toString());
				mOut.print("\t");
				mOut.print(sc.path);
				mOut.print("\t");
				mOut.print(sc.count);
				mOut.println();
				}
			mOut.flush();
			return 0;
			}
		catch(final Throwable err) {
			getLogger().error(err);
			return -1;
			}
		finally {
			for(SAMCount sw:tag2sw.values()) {
				sw.sfw.close();
				}
			}
		}
	}
