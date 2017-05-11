package com.github.lindenb.jvarkit.util.swing;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

import javax.swing.AbstractAction;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListModel;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;


 @SuppressWarnings("serial")
public class MultipleInputChooser extends AbstractFileChooser
	{
	public static final String EDITABLE="central.files.editable";
	public static final String PREDICATE="central.files.predicate";
	private JTextField textField;
	private AbstractAction setAction;
	private AbstractAction resetAction;
	private AbstractAction saveListAsAction;
	private DefaultListModel<String> listModel;
	private JList<String> list;
	private Predicate<String> predicate;
	public  MultipleInputChooser()
		{
		final JPanel top = new JPanel(new FlowLayout(FlowLayout.LEADING));
		super.add(top,BorderLayout.NORTH);
		
		this.listModel = new DefaultListModel<>();
		this.list = new JList<String>(this.listModel);
		final JScrollPane scroll = new JScrollPane(this.list);
		super.add(scroll,BorderLayout.CENTER);
		this.list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		
		
		this.textField = new JTextField(20);
		top.add(this.textField);
		this.textField.setEditable(true);
		this.setAction = new AbstractAction("[+]")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				chooseFiles();
				
				}
			};
		this.setAction.putValue(AbstractAction.LONG_DESCRIPTION, "Set File...");
		top.add(new JButton(this.setAction));
		this.resetAction = new AbstractAction("[-]")
			{
			@Override
			public void actionPerformed(ActionEvent e)
				{
				if(!textField.getText().trim().isEmpty())
					{	
					textField.setText("");
					return;
					}
				int indexes[]=list.getSelectedIndices();
				for(int i=indexes.length-1;i>=0;--i)
					{
					listModel.remove(indexes[i]);
					}
				}
			};
		this.resetAction.putValue(AbstractAction.LONG_DESCRIPTION, "Delete File...");
		this.resetAction.setEnabled(false);
		top.add(new JButton(this.resetAction));
		
		this.saveListAsAction = new AbstractAction("[!]")
			{
			@Override
			public void actionPerformed(ActionEvent e) {
				saveListAs();
				}
			};
		this.saveListAsAction.putValue(AbstractAction.LONG_DESCRIPTION, "Save List As...");
		this.saveListAsAction.setEnabled(false);
		top.add(new JButton(this.saveListAsAction));
		
		this.list.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent e)
				{
				resetAction.setEnabled(!list.isSelectionEmpty() || !textField.getText().trim().isEmpty());
				}
			});
		
		this.textField.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(textField.getText().trim().isEmpty()) return;
				if(addString(textField.getText().trim()))
					{
					textField.setText("");
					}
			}
		});
		
		}
	
	public JTextField getTextField()
		{
		return textField;
		}
	
	public ListModel<String> getListModel() {
		return listModel;
		}
	
	public List<String> getAsList()
		{
		final List<String> L = new ArrayList<>(this.listModel.size());
		for(int i=0;i< this.listModel.size();++i)
			{
			L.add(this.listModel.getElementAt(i));
			}
		return L;
		}
	
	private void chooseFiles()
		{
		final JFileChooser fc = new JFileChooser(PreferredDirectory.get(getClass()));
		fc.setFileFilter(getFileFilter());
		fc.setMultiSelectionEnabled(true);
		if(fc.showOpenDialog(this)!=JFileChooser.APPROVE_OPTION) return;
		for(File f : fc.getSelectedFiles())
			{
			addString(f.getPath());
			}
		textField.setText("");
		}
	
	private boolean addString(String s)
		{
		if(getPredicate()!=null && !getPredicate().test(s)) return false;
		for(int i=0;i< listModel.size();++i)
			{
			if(listModel.elementAt(i).equals(s))
				{
				return false;
				}
			}
		listModel.addElement(s);
		return true;
		}
	
	private void saveListAs()
		{
		final JFileChooser fc = new JFileChooser(PreferredDirectory.get(getClass()));
		if( fc.showSaveDialog(this)!=JFileChooser.APPROVE_OPTION) return;
		final File f=fc.getSelectedFile();
		if(f.exists() && JOptionPane.showConfirmDialog(this, "File "+f.getName()+" exists. Overwite ?", "Overwite ?", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
			{
			return;
			}
		PreferredDirectory.update(getClass(),f.getParentFile());
		PrintWriter pw =null;
		try {
			pw = new PrintWriter(f);
			for(int i=0;i< listModel.size();++i)
				pw.println(listModel.elementAt(i));
			pw.flush();
		} catch (Exception e) {
			JOptionPane.showMessageDialog(this, "Error "+e.getMessage());
			}
		finally
			{
			if(pw!=null) pw.close();
			}
		}
	
	public void setPredicate(Predicate<String> predicate)
		{
		Predicate<String> old=this.predicate;
		this.predicate = predicate;
		firePropertyChange(PREDICATE, old, predicate);
		}
	public Predicate<String> getPredicate() {
		return predicate;
		}
	}
