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
package org.mantoQSAR.gui;

import java.io.File;

public class GuiModel {

    private static GuiModel instance = null;
    File mantoFolder;
    String projectName;
    String descriptorSetting;
    String observationSetting;

    protected GuiModel() {

        File mainFolder = new File(System.getProperty("user.home"), "/Documents");

        if (!mainFolder.exists()) {
            mainFolder = new File(System.getProperty("user.home"), "/Dokumente");
        }

        File pF = new File(mainFolder, "/mantoProjects");

        if (!pF.exists()) {
            pF = new File(mainFolder, "/mantoQSARProjects");

        }

        if (!pF.exists()) {
            pF = mainFolder;
        }

        this.mantoFolder = pF;

        File[] projects = this.mantoFolder.listFiles();

        for (File project : projects) {
            if (project.isDirectory() == true) {
                this.projectName = project.getName();
                break;
            }
        }

        this.descriptorSetting = "descriptorSetting.json";
        this.observationSetting = "observationSetting.json";
    }

    public static GuiModel getInstance() {
        if (instance == null) {
            instance = new GuiModel();
        }

        return instance;

    }

    public void setSavedProject(File projectF) {
        File f = new File(projectF.getAbsolutePath());
        this.projectName = f.getName();
        this.mantoFolder = f.getParentFile();
    }

    public File getMantoProjectFolder() {
        return this.mantoFolder;
    }

    public void setMantoProjectFolder(File f) {
        this.mantoFolder = f;

    }

    public File getSavedProjectFolder() {
        return new File(this.mantoFolder.getAbsolutePath() + "/" + this.projectName);
    }

    public File getSavedDescriptorSetFile() {
        return new File(this.mantoFolder.getAbsolutePath() + "/" + this.projectName + "/set/" + this.descriptorSetting);

    }

    public File getSavedObservationSetFile() {
        return new File(this.mantoFolder.getAbsolutePath() + "/" + this.projectName + "/set/" + this.observationSetting);
    }

    public File getSavedDescriptorListFile() {
        return new File(this.mantoFolder.getAbsolutePath() + "/" + this.projectName + "/descriptorList.json");

    }
}
