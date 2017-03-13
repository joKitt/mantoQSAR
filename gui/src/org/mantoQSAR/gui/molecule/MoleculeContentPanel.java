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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.TransferHandler;
import javax.swing.border.BevelBorder;
import javax.swing.border.TitledBorder;
import org.jdesktop.swingx.JXCollapsiblePane;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.ObservationSet;
import org.mantoQSAR.core.util.ColorStatic;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MoleculeContentPanel extends JPanel {

    Logger logger;

    Molecule molecule;
    ObservationSet observationSet;
    MoleculeController moleculeController;

    private MouseListener mListener;
    private boolean isActive = false;
    public Integer id = 0;

    JPanel controlsTop;
    JPanel controlsDetail;

    private JXCollapsiblePane cPane;
    private TitledBorder border;

    private JTextField dropField; // field to drop file 

    private final Dimension textDim = new Dimension(150, 26);

    public MoleculeContentPanel(ObservationSet obsSet, Molecule molecule, Integer id) {

        this.logger = LoggerFactory.getLogger(MoleculeContentPanel.class);

        this.molecule = molecule;
        this.observationSet = obsSet;
        if (id != null) {
            this.id = id;
        }
        moleculeController = MoleculeController.getInstance();

        controlsTop = new JPanel(new GridBagLayout());
        controlsDetail = new JPanel(new GridBagLayout());

        initMoleculeContentPanel();

    }

    private void initMoleculeContentPanel() {

        // set mouse listener to get click
        mListener = new MouseListener() {
            @Override
            public void mouseClicked(MouseEvent me) {

                MoleculeController.getInstance().setActivePosition(id);
                // setActive(true);
            }

            @Override
            public void mousePressed(MouseEvent me) {
            }

            @Override
            public void mouseReleased(MouseEvent me) {
            }

            @Override
            public void mouseEntered(MouseEvent me) {
            }

            @Override
            public void mouseExited(MouseEvent me) {
            }
        };

        border = new TitledBorder(BorderFactory.createLineBorder(Color.LIGHT_GRAY), "New Entry");

        this.setBorder(border);
        this.setLayout(new BorderLayout());
        this.setBackground(Color.WHITE);

        cPane = new JXCollapsiblePane();
        cPane.setBackground(Color.WHITE);
        cPane.setLayout(new BorderLayout());

        controlsTop.setBackground(Color.WHITE);
        controlsDetail.setBackground(Color.WHITE);

        this.add(controlsTop, BorderLayout.PAGE_START);

        cPane.add(controlsDetail, BorderLayout.PAGE_START);
        cPane.setCollapsed(true);

        this.add(cPane, BorderLayout.CENTER);
        this.addMouseListener(mListener);
        this.composeForm();
    }

    private void composeForm() {

        this.controlsTop.removeAll();
        this.controlsDetail.removeAll();

        if (this.isActive == true) {
            controlsTop.setBackground(Color.YELLOW);
            controlsDetail.setBackground(Color.YELLOW);

        } else {
            controlsTop.setBackground(Color.WHITE);
            controlsDetail.setBackground(Color.WHITE);
        }

        this.setTitle(Integer.toString(this.id + 1));

        JTextField nameField = new JTextField(this.observationSet.getName());
        nameField.setPreferredSize(this.textDim);
        nameField.setMaximumSize(this.textDim);
        nameField.setBorder(BorderFactory.createLineBorder(Color.LIGHT_GRAY));

        final JLabel labelToggle = new JLabel("show");
        labelToggle.setForeground(Color.BLUE);

        labelToggle.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (cPane.isCollapsed()) {
                    cPane.setCollapsed(false);
                    labelToggle.setText("hide");
                } else {

                    cPane.setCollapsed(true);
                    labelToggle.setText("show");
                }
            }
        });

        controlsTop.add(labelToggle, newHideShowConstraints());
        controlsTop.add(nameField, newTextFieldConstraints(1));

        JCheckBox chkIsActive = new JCheckBox();
        chkIsActive.setOpaque(false);

        JCheckBox chkIsPredict = new JCheckBox();
        chkIsPredict.setOpaque(false);

        chkIsActive.setSelected(this.observationSet.isActive());

        if (this.observationSet.isActive() == true) {
            chkIsPredict.setSelected(this.observationSet.isPredict());
        } else {
            chkIsPredict.setEnabled(false);
        }

        controlsTop.add(chkIsActive, newCheckboxConstraints1());
        controlsTop.add(chkIsPredict, newCheckboxConstraints2());

        JLabel phLabel = new JLabel("pH value ");

        JTextField phField = new JTextField();
        phField.setText(String.valueOf(this.observationSet.getCondition().pH));
        phField.setPreferredSize(textDim);

        controlsDetail.add(phLabel, newLabelConstraints());
        controlsDetail.add(phField, newTextFieldConstraints(2));

        JLabel icLabel = new JLabel("ionic strength [mM] ");

        JTextField icField = new JTextField();
        icField.setText(String.valueOf(this.observationSet.getCondition().getIonicStrength()));
        icField.setPreferredSize(textDim);

        controlsDetail.add(icLabel, newLabelConstraints());
        controlsDetail.add(icField, newTextFieldConstraints(2));

        JLabel noteLabel = new JLabel("notes ");

        JTextArea noteArea = new JTextArea();
        noteArea.setEditable(true);
        noteArea.setSize(200, 200);
        noteArea.setPreferredSize(new Dimension(200, 200));
        noteArea.setText(this.observationSet.getNote());
        noteArea.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED, Color.LIGHT_GRAY, Color.LIGHT_GRAY));

        controlsDetail.add(noteLabel, newLabelConstraints());
        controlsDetail.add(noteArea, newTextFieldConstraints(2));

        JLabel labelSave = new JLabel("save entry");
        labelSave.setForeground(Color.BLUE);
        labelSave.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                System.out.println(ColorStatic.PURPLE + "Save not yet implemented." + ColorStatic.RESET);
            }
        });
    }

    private void addFileImport() {

        JLabel fileLabel = new JLabel("file ");

        dropField = new JTextField();
        dropField.setPreferredSize(textDim);

        try {
            dropField.setText(this.observationSet.getFile());
        } catch (Exception c) {

            dropField.setText("no valid file ");
            border.setBorder(BorderFactory.createLineBorder(Color.RED));

            JLabel dropLabel = new JLabel("Drop new structure file here", SwingConstants.CENTER);
            dropLabel.setForeground(Color.LIGHT_GRAY);
            controlsDetail.add(dropLabel, newTextFieldLastConstraints());
        }

        controlsDetail.add(fileLabel, newLabelConstraints());
        controlsDetail.add(dropField, newTextFieldConstraints(1));

        controlsDetail.setTransferHandler(new TransferHandler() {
            @Override
            public boolean canImport(TransferHandler.TransferSupport info) {
                return info.isDataFlavorSupported(DataFlavor.javaFileListFlavor);
            }

            @Override
            public boolean importData(TransferHandler.TransferSupport info) {
                if (!info.isDrop()) {
                    return false;
                }

                Transferable transferable = info.getTransferable();
                String importFileName = null;
                try {
                    @SuppressWarnings("unchecked")
                    List<File> fileList = (List<File>) transferable.getTransferData(DataFlavor.javaFileListFlavor);
                    Iterator<File> iterator = fileList.iterator();
                    while (iterator.hasNext()) {
                        File f = iterator.next();
                        importFileName = f.getAbsolutePath();
                    }

                    observationSet.setFile(importFileName);
                    dropField.setText(importFileName);

                } catch (UnsupportedFlavorException | IOException e) {
                    logger.error(e.getMessage());
                    return false;
                }
                logger.info("Importing " + importFileName);
                return true;
            }
        });
    }

    public void setTitle(String title) {
        border.setTitle(title);

    }

    public boolean isActive() {
        return isActive;
    }

    public void setActive(boolean isActive) {
        this.isActive = isActive;
    }

    private GridBagConstraints newConstraints() {
        GridBagConstraints c = new GridBagConstraints();
        c.insets = new Insets(2, 2, 2, 2);
        return c;
    }

    private GridBagConstraints newCheckboxConstraints1() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.LINE_END;
        c.weightx = 0.0;
        c.gridwidth = 1;
        return c;
    }

    private GridBagConstraints newCheckboxConstraints2() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.LINE_END;
        c.weightx = 0.0;
        c.gridwidth = 2;
        return c;
    }

    private GridBagConstraints newLabelHeaderConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.BASELINE_LEADING;
        c.weightx = 0.0;
        c.gridwidth = GridBagConstraints.REMAINDER;
        return c;
    }

    private GridBagConstraints newHideShowConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.FIRST_LINE_END;
        c.weightx = 0.0;
        c.gridwidth = GridBagConstraints.REMAINDER;
        return c;
    }

    private GridBagConstraints newLabelConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.BASELINE_TRAILING;
        c.weightx = 0.0;
        return c;
    }

    private GridBagConstraints newTextFieldConstraints(int n) {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.BASELINE_LEADING;
        c.weightx = 1.0;
        if (n == 2) {
            c.fill = GridBagConstraints.HORIZONTAL;
            c.gridwidth = GridBagConstraints.REMAINDER;
        }
        return c;
    }

    private GridBagConstraints newTextFieldLastConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.BASELINE_LEADING;
        c.weightx = 1.0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridwidth = GridBagConstraints.REMAINDER;
        return c;
    }

    private GridBagConstraints newBtnConstraints() {
        GridBagConstraints c = newConstraints();
        c.anchor = GridBagConstraints.NORTHEAST;
        c.weightx = 1.0;
        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridwidth = GridBagConstraints.REMAINDER;
        return c;
    }
}
