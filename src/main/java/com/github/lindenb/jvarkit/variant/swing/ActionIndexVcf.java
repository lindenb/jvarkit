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
package com.github.lindenb.jvarkit.variant.swing;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import javax.swing.AbstractAction;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileFilter;

import com.github.lindenb.jvarkit.swing.PreferredDirectory;
import com.github.lindenb.jvarkit.swing.ThrowablePane;
import com.github.lindenb.jvarkit.util.log.Logger;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.FileExtensions;
import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.IndexFactory.IndexType;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.vcf.VCFCodec;

@SuppressWarnings("serial")
public class ActionIndexVcf extends AbstractAction {
	private static final Logger LOG = Logger.of(ActionIndexVcf.class);

	public ActionIndexVcf() {
		super("Index VCF...");
		this.putValue(AbstractAction.SHORT_DESCRIPTION, "Index VCF files.");
		this.putValue(AbstractAction.LONG_DESCRIPTION, "Index one or more vcf/vcf.gz file(s)");
		}
	
	
	public static void indexVcf(final Path vcf,boolean overwrite) throws IOException {
		final boolean is_vcf_gz = vcf.getFileName().toString().toLowerCase().endsWith(FileExtensions.COMPRESSED_VCF) ;

		if(!overwrite) {
			final Path idxPath = is_vcf_gz ?
					Tribble.tabixIndexPath(vcf):
					Tribble.indexPath(vcf)
					;
			if(Files.exists(idxPath)) return;
			}
		
		final Path dir=vcf.getParent();
		try {
			if(dir!=null) IOUtil.assertDirectoryIsWritable(dir);
		} catch(SAMException err) {
			throw new IOException(err);
		}
		final SAMSequenceDictionary dict = SAMSequenceDictionaryExtractor.extractDictionary(vcf);
		
		final VCFCodec codec = new VCFCodec();
		LOG.warn("building index for "+vcf+".. please wait...");
		final Index index;
		if(is_vcf_gz) { 
			index = IndexFactory.createTabixIndex(vcf, codec, dict);
			}
		else
			{
			index = IndexFactory.createIndex(vcf, codec,IndexType.INTERVAL_TREE, dict);
			}
		index.writeBasedOnFeaturePath(vcf);
		}
	
	@Override
	public void actionPerformed(final ActionEvent e) {
		final JFileChooser fc = new JFileChooser(PreferredDirectory.get(ActionIndexVcf.class));
		fc.setMultiSelectionEnabled(true);
		fc.setFileFilter(new FileFilter() {
			@Override
			public String getDescription() {
				return "VCF files";
			}
			
			@Override
			public boolean accept(final File f) {
				if(!f.canRead()) return false;
				final String suff = f.getName().toLowerCase();
				return suff.endsWith(FileExtensions.COMPRESSED_VCF) || suff.endsWith(FileExtensions.VCF);
			}
		});
		final Object src = e.getSource();
		final Component parent = (src!=null && src instanceof Component ? Component.class.cast(src):null);
		if(fc.showOpenDialog(parent)!=JFileChooser.APPROVE_OPTION) {
			return;
			}
		final File[] vcfs = fc.getSelectedFiles();
		if(vcfs==null || vcfs.length==0) return;
		PreferredDirectory.update(ActionIndexVcf.class, vcfs[0]);
		this.setEnabled(false);
		try {
			final Thread t = new Thread() {
				@Override
				public void run() {
					int n_ok=0;
					for(File vcf:vcfs) {
						try {
							indexVcf(vcf.toPath(),true);
							n_ok++;
						} catch(Throwable err) {
							LOG.error(err);
						}
					}
					final int final_n_ok= n_ok;
					SwingUtilities.invokeLater(()-> {
						JOptionPane.showMessageDialog(parent,String.valueOf(final_n_ok)+"/"+vcfs.length+" "+
								(vcfs.length==1?"was":"were") +" successfuly indexed.",
								"Done.",
								JOptionPane.INFORMATION_MESSAGE
								);
						ActionIndexVcf.this.setEnabled(true);
					});
				}
			};
		t.start();
		} catch(final Throwable err) {
			ThrowablePane.show(parent, err);
			this.setEnabled(true);
		} finally {
			
		}
	}
}
