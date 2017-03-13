/* This file is part of mantoQSAR.

 mantoQSAR - Quantitative structure-activity relationship descriptor 
 calculation and modeling for biomolecules.
			
 Copyright (C) 2016  Jörg Kittelmann


 mantoQSAR is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, 
 or any later version.

 mantoQSAR is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with mantoQSAR. If not, see <http://www.gnu.org/licenses/>.
 */
package org.mantoQSAR.gui.molecule;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.FileFilter;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.border.EmptyBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.filechooser.FileSystemView;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import org.mantoQSAR.core.util.ColorStatic;
import org.mantoQSAR.gui.GuiModel;

public class ProjectTree extends JPanel {

    private FileSystemView fileSystemView;
    private final JTree tree;
    private DefaultTreeModel treeModel;

    private GuiModel guiModel;
    private MoleculeController moleculeController;

    public ProjectTree() {

        guiModel = GuiModel.getInstance();
        moleculeController = MoleculeController.getInstance();

        fileSystemView = FileSystemView.getFileSystemView();

        this.setLayout(new BorderLayout());
        DefaultMutableTreeNode root = new DefaultMutableTreeNode();
        treeModel = new DefaultTreeModel(root);
        TreeSelectionListener treeSelectionListener = new TreeSelectionListener() {

            @Override
            public void valueChanged(TreeSelectionEvent tse) {
                DefaultMutableTreeNode node
                        = (DefaultMutableTreeNode) tse.getPath().getLastPathComponent();
                File file = (File) node.getUserObject();

                if (!file.isFile()) {
                    showChildren(node, 2);
                }

                if (file.isFile()) {
                    showFileContent(node);
                }
            }

        };

        FileFilter directoryFilter = new FileFilter() {

            @Override
            public boolean accept(File file) {
                return file.isDirectory();
            }
        };

        File projectFolder = guiModel.getMantoProjectFolder();

        System.out.println(ColorStatic.PURPLE + "Project folder identified as " + projectFolder.toString() + ColorStatic.RESET);

        File[] projectRoot = projectFolder.listFiles(directoryFilter);

        DefaultMutableTreeNode node = new DefaultMutableTreeNode(projectFolder);
        root.add(node);
        showChildren(node, 3);

        MouseListener ml = new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                int selRow = tree.getRowForLocation(e.getX(), e.getY());
                TreePath selPath = tree.getPathForLocation(e.getX(), e.getY());
                if (selRow != -1) {
                    if (e.getClickCount() == 1) {

                    } else if (e.getClickCount() == 2) {
                        DefaultMutableTreeNode tn = (DefaultMutableTreeNode) tree.getSelectionPath().getLastPathComponent();
                        File f = (File) tn.getUserObject();
                        moleculeController.openProject(f);

                    }
                }
            }
        };

        tree = new JTree(treeModel);
        tree.setRootVisible(false);
        tree.addTreeSelectionListener(treeSelectionListener);

        tree.addMouseListener(ml);

        tree.setCellRenderer(new FileTreeCellRenderer());

        tree.expandRow(0);
        tree.setBorder(new EmptyBorder(5, 5, 5, 5));

        JScrollPane treeScroll = new JScrollPane(tree);
        add(treeScroll, BorderLayout.CENTER);
    }

    public void showRootFile() {
        tree.setSelectionInterval(0, 0);
    }

    private TreePath findTreePath(File find) {
        for (int ii = 0; ii < tree.getRowCount(); ii++) {
            TreePath treePath = tree.getPathForRow(ii);
            Object object = treePath.getLastPathComponent();
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) object;
            File nodeFile = (File) node.getUserObject();

            if (nodeFile == find) {
                return treePath;
            }
        }
        return null;
    }

    private void showErrorMessage(String errorMessage, String errorTitle) {
        JOptionPane.showMessageDialog(
                this,
                errorMessage,
                errorTitle,
                JOptionPane.ERROR_MESSAGE
        );
    }

    private void showChildren(final DefaultMutableTreeNode node, final Integer depth) {

        if (depth <= 0) {
            return;
        }
        File file = (File) node.getUserObject();
        if (file.isDirectory()) {
            File[] files = fileSystemView.getFiles(file, true); //!!

            node.removeAllChildren();

            for (File child : files) {
                DefaultMutableTreeNode childNode = new DefaultMutableTreeNode(child);

                if (!child.isDirectory() & depth > 2) {
                    // do not show files in main directory
                } else {
                    node.add(childNode);
                }

                if (child.isDirectory()) {
                    showChildren(childNode, depth - 1);
                }
            }
        }
    }

    private void showFileContent(DefaultMutableTreeNode node) {

        File file = (File) node.getUserObject();
        String extension = "";
        String fileName = file.toString();

        int i = fileName.lastIndexOf('.');
        if (i > 0) {
            extension = fileName.substring(i + 1);
        }

        if ("pqr".equals(extension) || "pdb".equals(extension)) {

        }

    }

    class FileTreeCellRenderer extends DefaultTreeCellRenderer {

        private final FileSystemView fileSystemView;
        private final JLabel label;

        private final ImageIcon folderCollapseImage = new ImageIcon(this.getClass().getResource("/org/mantoQSAR/gui/resources/icons/folder-open-img16x16.png"));
        private final ImageIcon folderExpandImage = new ImageIcon(this.getClass().getResource("/org/mantoQSAR/gui/resources/icons/folder-img16x16.png"));
        private final ImageIcon fileImage = new ImageIcon(this.getClass().getResource("/org/mantoQSAR/gui/resources/icons/text-x-generic-img16x16.png"));

        public FileTreeCellRenderer() {
            label = new JLabel();
            label.setOpaque(true);
            fileSystemView = FileSystemView.getFileSystemView();
        }

        @Override
        public Component getTreeCellRendererComponent(
                JTree tree,
                Object value,
                boolean selected,
                boolean expanded,
                boolean leaf,
                int row,
                boolean hasFocus) {

            DefaultMutableTreeNode node = (DefaultMutableTreeNode) value;
            File file = (File) node.getUserObject();

            label.setIcon(fileImage);

            if (file.isDirectory()) {
                label.setIcon(folderCollapseImage);
            }
            if (expanded) {
                label.setIcon(folderExpandImage);
            }
            if (file.isDirectory()) {

                label.setText(fileSystemView.getSystemDisplayName(file));
            } else {

                label.setText(file.getName());

            }
            label.setToolTipText(file.getPath());

            if (selected) {
                label.setBackground(backgroundSelectionColor);
                label.setForeground(textSelectionColor);
            } else {
                label.setBackground(backgroundNonSelectionColor);
                label.setForeground(textNonSelectionColor);
            }

            return label;
        }

    }

}
