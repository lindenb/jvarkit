package com.github.lindenb.jvarkit.util.swing;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

import javax.swing.AbstractAction;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileFilter;

@SuppressWarnings("serial")
public class FilePeeker extends JPanel
	{
	private File file=null;
	private FileFilter fileFilter=null;
	private String name="File";
	private JLabel label;
	private JTextField tField=null;
	private JList<File> fileList=null;
	private boolean multiple;
	
	public FilePeeker()
		{
		this(false);
		}
	
	public FilePeeker(boolean multiple)
		{
		super(new BorderLayout());
		this.multiple=multiple;
		this.label=new JLabel(name,JLabel.TRAILING);

		AbstractAction action=new AbstractAction("Choose") {
			@Override
			public void actionPerformed(ActionEvent e) {
				choose();
				}
			};
		this.getActionMap().put("choose.file", action);
		
		if(multiple)
			{
			JPanel top=new JPanel(new FlowLayout());
			add(top,BorderLayout.NORTH);
			top.add(label);
			
			this.fileList=new JList<File>(new DefaultListModel<File>());
			add(new JScrollPane(this.fileList),BorderLayout.CENTER);

			
			action=new AbstractAction("Remove")
				{
				@Override
				public void actionPerformed(ActionEvent e) {
					DefaultListModel<File> m=(DefaultListModel<File>)fileList.getModel();
					for(File f:fileList.getSelectedValuesList())
						{
						m.removeElement(f);
						}
					}
				};
			this.getActionMap().put("remove.file", action);
			action.setEnabled(false);
			this.fileList.getSelectionModel().addListSelectionListener(new ListSelectionListener()
				{
				@Override
				public void valueChanged(ListSelectionEvent e) {
					getActionMap().get("remove.file").setEnabled(!fileList.isSelectionEmpty());
					}
				});
			top.add(new JButton(this.getActionMap().get("choose.file")));
			top.add(new JButton(this.getActionMap().get("remove.file")));
			}
		else
			{
			add(label,BorderLayout.WEST);
			this.tField=new JTextField(20);
			this.tField.setEditable(false);
			add(this.tField,BorderLayout.CENTER);
			JButton but=new JButton(this.getActionMap().get("choose.file"));
			but.setText("Choose....");
			add(but,BorderLayout.EAST);
			}
		}
	
	public void setFiles(List<File> files)
		{
		if(multiple)
			{
			DefaultListModel<File> m=(DefaultListModel<File>)fileList.getModel();
			m.clear();
			for(File f:files)
				{
				if(!m.contains(f)) m.addElement(f);
				}
			}
		else
			{
			setFile(files.isEmpty()?null:files.get(0));
			}
		}
	
	public void setFile(File f)
		{
		if(multiple)
			{
			if(f==null)
				{
				setFiles(new ArrayList<File>());
				}
			else
				{
				setFiles(Collections.singletonList(f));
				}
			}
		else
			{
			this.tField.setText(this.file.getPath());
			this.tField.setCaretPosition(0);
			}
		}
	
	public void setFileFilter(FileFilter fileFilter)
		{
		this.fileFilter = fileFilter;
		}
	
	public FileFilter getFileFilter()
		{
		return fileFilter;
		}
	
	public FilePeeker filter(FileFilter filter)
		{
		setFileFilter(filter);
		return this;
		}
	
	public void choose()
		{
		Preferences prefs=Preferences.userNodeForPackage(getClass());
		String dirStr=prefs.get("last.directory."+name, null);
		File lastDir=null;
		if(dirStr!=null) lastDir=new File(dirStr);
		JFileChooser chooser=new JFileChooser(lastDir);
		if(fileFilter!=null) chooser.setFileFilter(fileFilter);
		chooser.setMultiSelectionEnabled(this.multiple);
		if(chooser.showOpenDialog(this)!=JFileChooser.APPROVE_OPTION)
			{
			return;
			}
		if(multiple)
			{
			DefaultListModel<File> m=(DefaultListModel<File>)fileList.getModel();

			for(File f:chooser.getSelectedFiles())
				{
				if(m.indexOf(f)!=-1) continue;
				m.addElement(f);
				}
			}
		else
			{
			File f=chooser.getSelectedFile();
			setFile(f);
			if(this.file.getParentFile()!=null)
				{
				prefs.put("last.directory."+name,
						this.file.getParentFile().getPath());
				try { prefs.sync();}catch(BackingStoreException err){}
				}
			}
		}
	public File getFile()
		{
		if(multiple)
			{
			DefaultListModel<File> m=(DefaultListModel<File>)fileList.getModel();
			if(m.isEmpty()) return null;
			return m.get(0);
			}
		else
			{
			return this.file;
			}
		}
	
	/** return a copy of the selected files */
	public List<File> getFiles()
		{
		if(multiple)
			{
			DefaultListModel<File> m=(DefaultListModel<File>)fileList.getModel();
			ArrayList<File> array=new ArrayList<File>(m.getSize());
			for(int i=0;i< m.getSize();++i)
				{
				array.add(m.get(i));
				}
			return array;
			}
		else
			{
			if( this.file==null) return Collections.emptyList();
			return Collections.singletonList(this.file);
			}
		}
	
	}
