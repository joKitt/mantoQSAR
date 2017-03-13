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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.IOException;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
import org.jdesktop.swingx.JXTable;
import org.mantoQSAR.core.ImportExport;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.util.ColorStatic;
import org.mantoQSAR.gui.Statics;
import org.mantoQSAR.gui.GuiModel;
import org.openide.windows.TopComponent;
import org.openide.windows.WindowManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MoleculeController implements ActionListener {

    final private PropertyChangeSupport changes;

    private static MoleculeController instance = null;
    Screen screen;
    Logger logger;
    GuiModel guiModel;

    JXTable descriptorTable;
    private int selectObs = 0;
    private int selectDescr = 30;
    private Integer orientationID = null;

    JmolPanel viewPanel = null;

    protected MoleculeController() {

        this.changes = new PropertyChangeSupport(this);
        this.logger = LoggerFactory.getLogger(MoleculeController.class);
        this.screen = Screen.getInstance();
        this.guiModel = GuiModel.getInstance();
        this.screen.setProjectPath(this.guiModel.getSavedProjectFolder().toString());
        this.descriptorTable = null;

    }

    public static MoleculeController getInstance() {

        if (instance == null) {
            instance = new MoleculeController();
        }
        return instance;
    }

    @Override
    public void actionPerformed(ActionEvent e) {

        logger.debug("Action call: " + e.getActionCommand());

        if (e.getActionCommand().equalsIgnoreCase(Statics.ACTION_CALC_DESCRIPTOR)) {
            screen.calcDescriptorAll();
        }

        if (e.getActionCommand().equalsIgnoreCase(Statics.ACTION_LOAD_DESCRIPTORSET)) {

            screen.loadDescriptorSetting(guiModel.getSavedDescriptorSetFile().toString());
        }

        if (e.getActionCommand().equalsIgnoreCase(MoleculeGuiStatic.LOAD_DESCRIPTOR_VALUES)) {
            logger.info("loading descriptor values from file " + guiModel.getSavedDescriptorListFile().toString());
            screen.loadDescriptorList(guiModel.getSavedDescriptorListFile().toString());
        }

        if (e.getActionCommand().equalsIgnoreCase(Statics.ACTION_LOAD_PROJECT)) {

            JFileChooser fileChooser = new JFileChooser(guiModel.getSavedProjectFolder());
            File selectedFolder = fileChooser.getSelectedFile();
            guiModel.setSavedProject(selectedFolder);

            logger.info("loading observation settings from file " + guiModel.getSavedObservationSetFile().toString());
            screen.loadObservationSetting(guiModel.getSavedObservationSetFile().toString());

            logger.info("loading descriptor settings from file " + guiModel.getSavedDescriptorSetFile().toString());
            screen.loadDescriptorSetting(guiModel.getSavedDescriptorSetFile().toString());

            logger.info("loading descriptor values from file " + guiModel.getSavedDescriptorListFile().toString());
            screen.loadDescriptorList(guiModel.getSavedDescriptorListFile().toString());

        }

        if (e.getActionCommand().equalsIgnoreCase(Statics.ACTION_SHOW_SURFACE)) {
            System.out.println("show surface fired in MoleculeController");

            this.setActiveDescriptor(21);
        }

        if (e.getActionCommand().equalsIgnoreCase(MoleculeGuiStatic.SAVE_DESCRIPTOR_TABLE)) {
            System.out.println("save descriptor table");
            this.saveDescriptor(this.descriptorTable);
        }

    }

    public Screen getScreen() {
        return screen;
    }

    void saveDescriptor(JXTable table) {

        ImportExport io = new ImportExport();
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.addChoosableFileFilter(new FileNameExtensionFilter("Microsoft Excel 97/2000/XP/2003 (.xls)", "xls"));
        fileChooser.addChoosableFileFilter(new FileNameExtensionFilter("Microsoft Excel 2007/2010/2013 XML (.xlsx)", "xlsx"));
        fileChooser.addChoosableFileFilter(new FileNameExtensionFilter("Text CSV (.csv)", "csv"));
        fileChooser.setAcceptAllFileFilterUsed(false);
        fileChooser.setSelectedFile(new File("Untitled 1"));

        TopComponent tc = WindowManager.getDefault().findTopComponent("DescriptorValueTopComponent");
        if (tc == null) {
            tc = new DescriptorValueTopComponent();
            tc.open();
        }

        int status = fileChooser.showSaveDialog(tc);

        if (status == JFileChooser.APPROVE_OPTION) {

            String[] ext = ((FileNameExtensionFilter) fileChooser.getFileFilter()).getExtensions();
            File file = fileChooser.getSelectedFile();
            String nameLower = file.getName().toLowerCase();

            for (String e : ext) {
                if (!nameLower.endsWith("." + e.toLowerCase())) {
                    file = new File(file.toString() + "." + ext[0]);
                }
            }

            System.out.println(ColorStatic.GREEN + "save table to " + file.toString() + ColorStatic.RESET);

            if (ext[0].equalsIgnoreCase("xls") || ext[0].equalsIgnoreCase("xlsx")) {
                try {
                    io.saveTableToExcel(table, file);
                } catch (IOException ex) {
                    System.out.println(ColorStatic.RED + ex.getMessage() + ColorStatic.RESET);
                }
            }

            if (ext[0].equalsIgnoreCase("csv")) {
                try {
                    io.saveTableToCSV(table, file);
                } catch (IOException ex) {
                    System.out.println(ColorStatic.RED + ex.getMessage() + ColorStatic.RESET);
                }
            }
        }
    }

    void openProjectDialog() {
        JFileChooser fileChooser = new JFileChooser(guiModel.getMantoProjectFolder());

        TopComponent tc = WindowManager.getDefault().findTopComponent("MoleculeSettingTopComponent");
        if (tc == null) {
            tc = new MoleculeSettingTopComponent();
            tc.open();
        }

        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        int status = fileChooser.showOpenDialog(tc);
        if (status == JFileChooser.APPROVE_OPTION) {
            File f = fileChooser.getSelectedFile();
            openProject(f);
        }

    }

    @Deprecated
    void openProject() {
        JFileChooser fileChooser = new JFileChooser(guiModel.getMantoProjectFolder());

        TopComponent tc = WindowManager.getDefault().findTopComponent("MoleculeSettingTopComponent");
        if (tc == null) {
            tc = new MoleculeSettingTopComponent();
            tc.open();
        }

        fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        int status = fileChooser.showOpenDialog(tc);

        if (status == JFileChooser.APPROVE_OPTION) {
            File f = fileChooser.getSelectedFile();
            if (!f.isDirectory()) {
                f = f.getParentFile();

                System.out.println(ColorStatic.PURPLE + f.toString() + ColorStatic.RESET);

                String folder = f.toString();
                if (folder.contains("model")
                        || folder.contains("set")
                        || folder.contains("structure")) {
                    f = f.getParentFile();
                }
            }

            guiModel.setSavedProject(f);
            screen.setProjectPath(f.toString());

        }

        this.screen.loadObservationSetting(guiModel.getSavedObservationSetFile().toString());
        this.screen.loadDescriptorSetting(guiModel.getSavedDescriptorSetFile().toString());
    }

    void openProject(File f) {

        if (!f.isDirectory()) {
            f = f.getParentFile();

            System.out.println(ColorStatic.PURPLE + f.toString() + ColorStatic.RESET);

            String folder = f.toString();
            if (folder.contains("model")
                    || folder.contains("set")
                    || folder.contains("structure")) {
                f = f.getParentFile();
            }
        }

        guiModel.setSavedProject(f);
        screen.setProjectPath(f.toString());

        this.screen.loadObservationSetting(guiModel.getSavedObservationSetFile().toString());
        this.screen.loadDescriptorSetting(guiModel.getSavedDescriptorSetFile().toString());
        this.screen.loadDescriptorList(guiModel.getSavedDescriptorListFile().toString());
    }

    public Integer getActivePosition() {
        return this.selectObs;
    }

    public void setActivePosition(int i) {
        Integer oldObs = this.selectObs;
        this.selectObs = i;
        System.out.println(ColorStatic.PURPLE + "active molecule set to id " + i + ColorStatic.RESET);
        changes.firePropertyChange(MoleculeGuiStatic.CHANGE_ACTIVE_OBSERVATION, oldObs.intValue(), this.selectObs);
    }

    public Integer getActiveDescriptor() {
        return this.selectDescr;
    }

    public void setActiveDescriptor(int i) {
        Integer oldDescr = this.selectDescr;
        this.selectDescr = i;
        System.out.println(ColorStatic.PURPLE + "active descriptor set to id " + i + ColorStatic.RESET);
        changes.firePropertyChange(MoleculeGuiStatic.CHANGE_ACTIVE_DESCRIPTOR, oldDescr.intValue(), this.selectDescr);
    }

    public JmolPanel getViewPanel() {
        if (this.viewPanel == null) {
            this.viewPanel = new JmolPanel();
        }
        return this.viewPanel;
    }

    public Integer getOrientationID() {
        return orientationID;
    }

    public void setOrientationID(Integer ID) {
        Integer oldOrientationID = this.orientationID;
        this.orientationID = ID;
        changes.firePropertyChange(MoleculeGuiStatic.CHANGE_ACTIVE_ORIENTATION, oldOrientationID, this.orientationID);
    }

    public void addPropertyChangeListener(PropertyChangeListener listener) {
        changes.addPropertyChangeListener(listener);
    }

    public void removePropertyChangeListener(PropertyChangeListener listener) {
        changes.removePropertyChangeListener(listener);
    }
}
