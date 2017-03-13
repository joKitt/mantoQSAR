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
package org.mantoQSAR.gui.model;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JPanel;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.mantoQSAR.core.ImportExport;
import org.mantoQSAR.core.model.ModelGroup;
import org.mantoQSAR.core.model.ModelList;
import org.mantoQSAR.core.model.ModelStatic;
import org.mantoQSAR.core.util.ColorStatic;
import org.mantoQSAR.gui.molecule.DescriptorValueTopComponent;
import org.mantoQSAR.gui.resources.CollapsiblePane;
import org.openide.windows.TopComponent;
import org.openide.windows.WindowManager;
import org.slf4j.LoggerFactory;
import org.slf4j.Logger;

public class ModelResultGraphSetting extends JPanel {

    final private Logger logger;
    private final ModelController modelControl;
    private ModelGroup modelGroup;
    private ModelList modelList;

    private final ModelStatisticPane mStatPane;

    private final ModelResultGraph modelGraph;

    public ModelResultGraphSetting(ModelResultGraph mG) {

        logger = LoggerFactory.getLogger(ModelResultGraphSetting.class);

        modelControl = ModelController.getInstance();
        modelGroup = modelControl.getModelGroup();
        modelList = modelGroup.getModelList();

        this.modelGraph = mG;

        JPanel groupPanel = new JPanel();
        groupPanel.setLayout(new BoxLayout(groupPanel, BoxLayout.Y_AXIS));

        this.setLayout(new BorderLayout());

        groupPanel.setBackground(Color.WHITE);

        mStatPane = new ModelStatisticPane();
        groupPanel.add(mStatPane);

        groupPanel.add(getSettingPane());
        groupPanel.add(new JPanel());
        this.add(groupPanel, BorderLayout.CENTER);

        initListeners();
    }

    protected final void initListeners() {

        modelControl.addPropertyChangeListener(new PropertyChangeListener() {

            @Override
            public void propertyChange(PropertyChangeEvent evt) {

                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CALCULATED)) {
                    updateStatistics();
                }
                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CHANGED)) {
                    updateStatistics();
                }
            }
        });

        modelList.addPropertyChangeListener(new PropertyChangeListener() {

            @Override
            public void propertyChange(PropertyChangeEvent evt) {

                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CALCULATED)) {
                    updateStatistics();
                }
                if (evt.getPropertyName().equalsIgnoreCase(ModelStatic.EVENT_MODEL_CHANGED)) {
                    updateStatistics();
                }
            }
        });
    }

    private CollapsiblePane getSettingPane() {

        CollapsiblePane cSettingPane = new CollapsiblePane();
        cSettingPane.setCollapsed(false);
        cSettingPane.setTitle("model setting");

        JPanel content = cSettingPane.getMainPanel();

        // check box to toggle detail appearance
        JCheckBox chkShowDetail = new JCheckBox("model detail");
        chkShowDetail.setBackground(Color.WHITE);
        chkShowDetail.setSelected(false);
        chkShowDetail.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {

                Boolean selected = ((JCheckBox) e.getSource()).isSelected();
                modelGraph.setModelDetail(selected);
            }
        });

        content.add(chkShowDetail);

        JCheckBox chkLogScale = new JCheckBox("log scale");
        chkLogScale.setBackground(Color.WHITE);
        chkLogScale.setSelected(false);
        chkLogScale.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {

                System.out.println(ColorStatic.PURPLE + "toggled model log scale");
                Boolean selected = ((JCheckBox) e.getSource()).isSelected();
                modelGraph.setLogScale(selected);
            }
        });
        content.add(chkLogScale);

        JButton btnSave = new JButton("save figure");
        btnSave.setBackground(Color.WHITE);
        btnSave.addActionListener(new ActionListener() {

            @Override
            public void actionPerformed(ActionEvent e) {

                ImportExport io = new ImportExport();
                JFileChooser fileChooser = new JFileChooser();
                fileChooser.addChoosableFileFilter(new FileNameExtensionFilter("Portable Network Graphic (.png)", "png"));
                fileChooser.setAcceptAllFileFilterUsed(false);
                fileChooser.setSelectedFile(new File("Untitled 1"));
                TopComponent tc = WindowManager.getDefault().findTopComponent("DescriptorValueTopComponent");
                if (tc == null) {
                    tc = new DescriptorValueTopComponent();
                    tc.open();
                }

                int status = fileChooser.showSaveDialog(tc);
                if (status == JFileChooser.APPROVE_OPTION) {
                    String[] extList = ((FileNameExtensionFilter) fileChooser.getFileFilter()).getExtensions();
                    File file = fileChooser.getSelectedFile();

                    String nameLower = file.getName().toLowerCase();

                    for (String ext : extList) {
                        if (!nameLower.endsWith("." + ext.toLowerCase())) {
                            file = new File(file.toString() + "." + extList[0]);
                        }
                    }

                    System.out.println(ColorStatic.GREEN + "save model figure to " + file.toString() + ColorStatic.RESET);
                    io.savePanelToImage(modelGraph, file);

                    modelGroup.displayModelData();

                }
            }
        });

        cSettingPane.getTopPane().add(btnSave);
        return cSettingPane;
    }

    public void updateStatistics() {

        mStatPane.setRSquare(modelGroup.getRSquare());
        mStatPane.setPredRSquare(modelGroup.getPredRSquare());
        this.revalidate();
        this.repaint();
    }

}
