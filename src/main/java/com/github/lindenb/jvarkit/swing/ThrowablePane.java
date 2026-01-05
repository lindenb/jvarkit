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
package com.github.lindenb.jvarkit.swing;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.HeadlessException;
import java.awt.Toolkit;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import javax.swing.BorderFactory;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableModel;


/**
 * Panel used to display the content of an exception
 * @author pierre
 *
 */
@SuppressWarnings("serial")
public class ThrowablePane extends JPanel
	{

    /**
     * Constructor with a throwable
     * @param throwable
     */
    public ThrowablePane(Throwable throwable,String message)
        {
        super(new BorderLayout());
        if(message==null) message=throwable.getLocalizedMessage();
        if(message==null) message=throwable.getMessage();
        Dimension screen= Toolkit.getDefaultToolkit().getScreenSize();
        setPreferredSize(new Dimension(screen.width/2,screen.height/2));
        setBorder(BorderFactory.createLineBorder(Color.RED,2));
        List<StackTraceElement> rows = new ArrayList<>();
        for(final StackTraceElement ste: throwable.getStackTrace())
        	{
        	rows.add(ste);
        	}
        final TableModel tableModel = new ColumnDefTableModel<StackTraceElement>(
        		Arrays.asList(
        		new ColumnDef<StackTraceElement>("Class",String.class,E->E.getClassName()),
        		new ColumnDef<StackTraceElement>("Method",String.class,E->E.getMethodName()),
        		new ColumnDef<StackTraceElement>("File",String.class,E->E.getFileName()),
        		new ColumnDef<StackTraceElement>("Line",Integer.class,E->E.getLineNumber())
        		),rows);
        
        JTabbedPane tabbedPane=new JTabbedPane(JTabbedPane.SCROLL_TAB_LAYOUT);
        this.add(tabbedPane,BorderLayout.CENTER);
        JTextField tf= new JTextField(throwable.getClass()+" : "+message);
        tf.setToolTipText(tf.getText());
        tf.setEditable(false);
        tf.setCaretPosition(0);
        tf.setForeground(Color.RED);
        tf.setFont(new Font("Dialog",Font.PLAIN,24));
        tabbedPane.addTab("Message", tf);
         
        JPanel pane= new JPanel(new GridLayout(0,1));
        tabbedPane.addTab("Trace", pane);
        StringWriter statckTrace= new StringWriter();
        throwable.printStackTrace(new PrintWriter(statckTrace));
        JTextArea msg= new JTextArea(throwable.getClass().getName()+
                "\n"+
                statckTrace.toString());
        msg.setEditable(false);
        msg.setBackground(Color.WHITE);
        msg.setLineWrap(true);
        JScrollPane scroll=new JScrollPane(msg);
        scroll.setBorder(BorderFactory.createTitledBorder("An exception of type "+
                throwable.getClass().getName() +
                " has occured"));
        scroll.setPreferredSize(new Dimension(200,200));
        pane.add(scroll);

        JTable table = new JTable(tableModel);
        DefaultTableCellRenderer render= new DefaultTableCellRenderer()
        	{
			private static final long serialVersionUID = 1L;

			@Override
        	public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
				{
        		Component c= super.getTableCellRendererComponent(table, value, isSelected, hasFocus,
        				row, column);
        		Object pack= table.getModel().getValueAt(row, 0);
        		
        		if(!isSelected &&
        			pack!=null && 
        			!pack.toString().contains("com.github.lindenb") )
        			{
        			setBackground(Color.GREEN);
        			}
        		else if(isSelected)
        			{
        			setBackground(Color.RED);
        			}
        		else
	    			{
	    			setBackground(Color.WHITE);
	    			}
        		return c;
        		}
        	};
        render.setOpaque(true);
        for(int i=0;i< table.getColumnModel().getColumnCount();++i)
        	{
        	table.getColumnModel().getColumn(i).setCellRenderer(render);
        	}
        table.setBackground(Color.WHITE);
        scroll=new JScrollPane(table);
        scroll.setBorder(BorderFactory.createTitledBorder("Stack Trace"));
        pane.add(scroll);
        
        /** try to display the file */
        for(StackTraceElement ste: rows)
        	{
        	String classname = ste.getClassName();
        	if(classname==null) continue;
        	classname.replace('.', '/');
        	classname="/"+classname+".java";
        	InputStream in = null;
        	try {
        		in= getClass().getResourceAsStream(classname);
            	if(in==null) continue;
            	String javaFile = IOUtil.readFully(in);
				JTextArea area= new JTextArea(javaFile.toString());
				JPanel srcPane= new JPanel(new BorderLayout());
				srcPane.setBorder(BorderFactory.createTitledBorder("Source"));
				srcPane.add(new JScrollPane(area));
				tabbedPane.addTab("Source", srcPane);
				}
        	catch (Throwable e)
				{
				break;
				}
        	finally
        		{
        		CloserUtil.close(in);
        		}		
        	break;
        	}
        
        setBackground(Color.RED);
        pane.setOpaque(true);
        pane.setBackground(Color.RED);
        }

    /** display an alert via OptionPane.showMessageDialog */
    static public void show(Component owner, String message,Throwable t)
        {
        try {
			JOptionPane.showMessageDialog(owner,new ThrowablePane(t,message),"Error",JOptionPane.ERROR_MESSAGE,null);
			}
        catch (HeadlessException e)
			{
			t.printStackTrace(System.err);
			}
        }
    
    /** display an alert via OptionPane.showMessageDialog */
    static public void show(Component owner, Throwable t)
        {
    	show(owner,null,t);
        }
    

}