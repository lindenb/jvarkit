package com.github.lindenb.jvarkit.tools.misc;

import htsjdk.samtools.util.CloserUtil;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.PrintWriter;
import java.net.URLEncoder;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.AbstractAction;
import javax.swing.JButton;
import javax.swing.JDesktopPane;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.event.InternalFrameAdapter;
import javax.swing.event.InternalFrameEvent;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.xml.XMLConstants;
import javax.xml.namespace.QName;
import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;
import javax.xml.transform.stream.StreamSource;

import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.jcommander.Program;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.ns.RDF;

/*
## Example

```bash
$  java -jar dist/buildwpontology.jar
```

*/
@SuppressWarnings("serial")
@Program(name="buildwpontology",
	description="Build a simple RDFS/XML ontology from the Wikipedia Categories",
	keywords={"wikipedia","ontology","rdf","gui"}
	)
public class BuildWikipediaOntology extends Launcher
	{
	final static String RDFS="http://www.w3.org/2000/01/rdf-schema#";
	private static final Logger LOG = Logger.build(BuildWikipediaOntology.class).make();

	private static long ID_GENERATOR=0L;
	private class Frame extends JFrame
		{
		private JDesktopPane desktopPane;
		
		
		private class IFrame extends JInternalFrame
			{
			private JTree jtree;
			private File saveAs;
			
			private class CatNode extends DefaultMutableTreeNode
				{
				private long id=(++ID_GENERATOR);
				private String catName;
				CatNode(String catName)
					{
					super(catName,true);
					this.catName=catName;
					}
				
				
				
				boolean isMainNode()
					{
					return getIdByTerm(this.catName)==this.id;
					}
				
				private void delete()
					{
					CatNode parent= (CatNode)CatNode.this.getParent();
					if(parent==null) return;
					CatNode.this.removeFromParent();
					DefaultTreeModel model=(DefaultTreeModel)jtree.getModel();
					model.nodeStructureChanged(parent);
					}
				
				private void expand()
					{
					Set<String> terms=getSubCategories(CatNode.this.catName);
					for(int i=0;i< CatNode.this.getChildCount();++i)
						{
						CatNode c=CatNode.class.cast( CatNode.this.getChildAt(i));
						terms.remove(c.catName);
						}
					
					for(String term:terms)
						{
						CatNode newnode=new CatNode(term);
						CatNode.this.add(newnode);
						}
					DefaultTreeModel model=(DefaultTreeModel)jtree.getModel();
					model.nodeStructureChanged(CatNode.this);
					}
				@Override
				public String toString()
					{
					String s=this.catName;
					int colon=s.indexOf(':');
					if(colon!=-1) s=s.substring(colon+1);
					if(!isMainNode()) s+="*";
					return s;
					}
				}
			
			public IFrame(String name)
				{
				super(name,true,true,true,true);
				this.setDefaultCloseOperation(JInternalFrame.DISPOSE_ON_CLOSE);
				JPanel contentPanel=new JPanel(new BorderLayout(5,5));
				this.setContentPane(contentPanel);
				contentPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
				contentPanel.add(new JScrollPane(this.jtree=new JTree(new DefaultTreeModel(new CatNode(name)))));
				this.addInternalFrameListener(new InternalFrameAdapter()
					{
					@Override
					public void internalFrameOpened(InternalFrameEvent e) {
						removeInternalFrameListener(this);
						CatNode.class.cast(jtree.getModel().getRoot()).expand();
						}
					});
				
				

				
				JToolBar toolBar=new JToolBar();
				contentPanel.add(toolBar,BorderLayout.NORTH);
				AbstractAction action=new AbstractAction("Expand")
					{
					@Override
						public void actionPerformed(ActionEvent e) {
						TreePath selPaths[]=jtree.getSelectionPaths();
						if(selPaths==null || selPaths.length==0)return;
						for(TreePath selPath:selPaths)
							{
							CatNode node=(CatNode)selPath.getLastPathComponent();
							node.expand();
							}
						}
					};
				toolBar.add(new JButton(action));
				action.setEnabled(false);
				getActionMap().put("action.expand", action);
				
				action=new AbstractAction("Delete")
					{
					@Override
						public void actionPerformed(ActionEvent e) {
						TreePath selPaths[]=jtree.getSelectionPaths();
						if(selPaths==null || selPaths.length==0)return;
						for(TreePath selPath:selPaths)
							{
							CatNode node=(CatNode)selPath.getLastPathComponent();
							node.delete();
							}
						}
					};
				toolBar.add(new JButton(action));
				action.setEnabled(false);
				getActionMap().put("action.delete", action);

				
				jtree.getSelectionModel().addTreeSelectionListener(new TreeSelectionListener()
					{
					@Override
					public void valueChanged(TreeSelectionEvent e)
						{	
						TreePath path=jtree.getSelectionPath();
						CatNode node=null;
						if(path!=null) node=CatNode.class.cast(path.getLastPathComponent());
						getActionMap().get("action.expand").setEnabled(node!=null && node.isMainNode());
						getActionMap().get("action.delete").setEnabled(node!=null && node.getParent()!=null);
						}
					});
				
				
				JMenuBar menuBar=new JMenuBar();
				setJMenuBar(menuBar);
				JMenu menu=new JMenu("File");
				menuBar.add(menu);
				menu.add(new AbstractAction("Save As...")
					{
					@Override
					public void actionPerformed(ActionEvent e) {
						doMenuSaveAs();
						}
					});

				menu.add(new AbstractAction("Save")
					{
					@Override
					public void actionPerformed(ActionEvent e) {
						doMenuSave(saveAs);
						}
					});
				
				
				menu.add(new AbstractAction("Close")
					{
					@Override
					public void actionPerformed(ActionEvent e) {
						IFrame.this.setVisible(false);
						IFrame.this.dispose();
						}
					});

				}
			String getURL(String catName)
				{
				return	"http://en.wikipedia.org/wiki/"+catName.replace(' ', '_');
				}
			
			private void visit(XMLStreamWriter w,CatNode root) throws XMLStreamException
				{
				if(root.isMainNode())
					{
					w.writeStartElement("rdfs", "Class", RDFS);
					w.writeAttribute("rdf", RDF.NS, "about", getURL(root.catName));
					int colon=root.catName.indexOf(':');
					
					w.writeStartElement("rdfs", "label", RDFS);
					w.writeCharacters(colon==-1?root.catName:root.catName.substring(colon+1));
					w.writeEndElement();
					
					for(String parent: getParentTerms(root.catName))
						{
						w.writeEmptyElement("rdfs", "subClassOf", RDFS);
						w.writeAttribute("rdf", RDF.NS, "resource",
								getURL(parent)
								);
						}
					
					w.writeEndElement();
					w.writeCharacters("\n");
					}
				for(int i=0;i< root.getChildCount();++i)
					{
					CatNode c=CatNode.class.cast(root.getChildAt(i));
					visit(w,c);
					}

				}
			private void doMenuSave(File f)
				{
				if(f==null)
					{
					doMenuSaveAs();
					return;
					}
				PrintWriter pw=null;
				XMLStreamWriter w=null;
				try
					{
					
					pw=new PrintWriter(f,"UTF-8");
					XMLOutputFactory xmlfactory= XMLOutputFactory.newInstance();
					w= xmlfactory.createXMLStreamWriter(pw);
					w.writeStartDocument("UTF-8","1.0");
					w.writeStartElement("rdf", "RDF", RDF.NS);
					w.writeAttribute("xmlns", XMLConstants.XML_NS_URI, "rdf",RDF.NS);
					w.writeAttribute("xmlns", XMLConstants.XML_NS_URI,  "rdfs", RDFS);
					w.writeCharacters("\n");
					visit(w,CatNode.class.cast(jtree.getModel().getRoot()));
					
					w.writeEndElement();
					w.writeEndDocument();
					this.saveAs=f;
					}
				catch(Exception err)
					{
					JOptionPane.showMessageDialog(this, "Error"+err.getMessage());
					}
				finally
					{
					CloserUtil.close(w);
					CloserUtil.close(pw);
					}
				}
			private void doMenuSaveAs()
				{
				JFileChooser chooser=new JFileChooser(saveAs==null?null:saveAs.getParentFile());
				if(chooser.showSaveDialog(this)!=JFileChooser.APPROVE_OPTION) return;
				File f=chooser.getSelectedFile();
				if(f.exists() && JOptionPane.showConfirmDialog(this, f.getPath()+" exists. Overwite ?", "Confirm", JOptionPane.OK_CANCEL_OPTION, JOptionPane.WARNING_MESSAGE, null)!=JOptionPane.OK_OPTION)
					{
					return;
					}
				doMenuSave(chooser.getSelectedFile());
				}
			
			private Set<String> getParentTerms(CatNode root,String term,Set<String> parents)
				{
			
				for(int i=0;i< root.getChildCount();++i)
					{
					CatNode c=CatNode.class.cast(root.getChildAt(i));
					if(c.catName.equals(term)) parents.add(root.catName);
					getParentTerms(c,term,parents);
					}
				return parents;
				}
			Set<String> getParentTerms(String term)
				{
				if(jtree==null || jtree.getModel()==null) return Collections.emptySet();
				return getParentTerms(CatNode.class.cast(jtree.getModel().getRoot()),term,new HashSet<String>());
				}

			
			private long getIdByTerm(CatNode root,String term,long curr)
				{
				if(root.catName.equals(term) && (curr==-1L || curr>root.id))
					{
					curr=root.id;
					}
				for(int i=0;i< root.getChildCount();++i)
					{
					CatNode c=CatNode.class.cast(root.getChildAt(i));
					long n=getIdByTerm(c,term,curr);
					if(curr==-1 || curr>n ) curr=n;
					}
				return curr;
				}

			
			
			long getIdByTerm(String term)
				{
				if(jtree==null || jtree.getModel()==null) return -1;
				return getIdByTerm(CatNode.class.cast(jtree.getModel().getRoot()),term,-1L);
				}
			
			
			}
		
		Frame()
			{
			super(getProgramName());
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			this.addWindowListener(new WindowAdapter()
				{
				@Override
				public void windowOpened(WindowEvent e) {
					doMenuNewOntology();
					}
				});
			this.desktopPane=new JDesktopPane();
			
			JPanel contentPanel=new JPanel(new BorderLayout());
			setContentPane(contentPanel);
			contentPanel.add(desktopPane,BorderLayout.CENTER);
			
			JMenuBar menuBar=new JMenuBar();
			setJMenuBar(menuBar);
			JMenu menu=new JMenu("File");
			menuBar.add(menu);
			menu.add(new AbstractAction("Quit")
				{
				@Override
				public void actionPerformed(ActionEvent e) {
					Frame.this.setVisible(false);
					Frame.this.dispose();
					}
				});
			menu=new JMenu("Tools");
			menuBar.add(menu);
			menu.add(new AbstractAction("New Ontology")
				{
				@Override
				public void actionPerformed(ActionEvent e) {
					doMenuNewOntology();
					}
				});
			}
		private void doMenuNewOntology()
			{
			try {
				String s=JOptionPane.showInputDialog(this, "Select Root Category","Category:Science");
				if(s==null) return;
				if(!s.startsWith("Category:"))
					{
					JOptionPane.showMessageDialog(this, "Root should start with 'Category:'");
					return;
					}
				Dimension d=desktopPane.getSize();
				IFrame iframe=new IFrame(s);
				iframe.setBounds(50, 50, d.width-100, d.height-100);
				desktopPane.add(iframe);
				iframe.setVisible(true);
				}
			catch (Exception err)
				{
				LOG.error(err);
				}	
			}	
		
		private Set<String> getSubCategories(String term)
			{
			Set<String> set=new HashSet<String>();
			String cmcontinue=null;
			XMLEventReader r=null;
			QName title=new QName("title");
			QName att_cmcontinue=new QName("cmcontinue");

			try
				{
				XMLInputFactory xif=XMLInputFactory.newFactory();
				do
				{
				
				String url="http://en.wikipedia.org/w/api.php?action=query&list=categorymembers&cmnamespace=14&cmlimit=100&format=xml&cmtitle="+
						URLEncoder.encode(term,"UTF-8")+
						(cmcontinue==null?"":"&cmcontinue="+cmcontinue)
						;
				LOG.info(url);
				cmcontinue=null;
				r=xif.createXMLEventReader(new StreamSource(url));
				while(r.hasNext())
					{
					XMLEvent evt=r.nextEvent();
					if(!evt.isStartElement()) continue;
					StartElement E=evt.asStartElement();
					Attribute att=null;
					if(E.getName().getLocalPart().equals("cm") && (att=E.getAttributeByName(title))!=null)
						{
						set.add(att.getValue());
						}
					else if(E.getName().getLocalPart().equals("categorymembers") && (att=E.getAttributeByName(att_cmcontinue))!=null)
						{
						cmcontinue=att.getValue();
						}
					}
				} while(cmcontinue!=null);
				}
			catch (Exception err)
				{
				LOG.error(err);
				}
			return set;
			}

		
		}
	
	private BuildWikipediaOntology()
		{
		
		}
	
	
	
	@Override
	public int doWork(List<String> args)
		{
		final Frame frame=new Frame();
		
		
		try
			{
			Dimension d=Toolkit.getDefaultToolkit().getScreenSize();
			frame.setBounds(100, 100, d.width-200, d.height-200);
			SwingUtilities.invokeAndWait(new Runnable()
				{
				@Override
				public void run() {
					frame.setVisible(true);
				}
			});
			return 0;
			}
		catch(Exception err)
			{
			LOG.error(err);
			return -1;
			}
		finally
			{
			
			}
		}
	public static void main(String[] args) {
		new BuildWikipediaOntology().instanceMain(args);
		}
	}
