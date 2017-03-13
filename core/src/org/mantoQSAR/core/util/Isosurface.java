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
package org.mantoQSAR.core.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import org.biojava.nbio.structure.Structure;
import java.io.File;
import java.io.RandomAccessFile;
import javax.swing.JPanel;
import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolStatusListener;
import org.jmol.api.JmolViewer;
import org.mantoQSAR.core.Screen;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Isosurface implements Runnable {

    Double resolution;
    Double sizeOfSphere;
    String typeOfSurface;
    Structure structure = null;
    private final long randID;
    String filepath;

    Thread t;

    public Isosurface() {
        this.resolution = 1.1;
        this.sizeOfSphere = 1.4;
        this.typeOfSurface = "sasurface";
        this.randID = generateRandom(12);
        this.filepath = Screen.getInstance().getProjectPath();

    }

    public static long generateRandom(int length) {
        Random random = new Random();
        char[] digits = new char[length];
        digits[0] = (char) (random.nextInt(9) + '1');
        for (int i = 1; i < length; i++) {
            digits[i] = (char) (random.nextInt(10) + '0');
        }
        return Long.parseLong(new String(digits));
    }

    public List<Double[]> getSurface(Structure struct) {

        File tempDirect = new File(filepath + "/temp");
        if (!tempDirect.exists()) {
            System.out.println("creating directory: " + tempDirect);

            try {
                tempDirect.mkdir();
                System.out.println("directory created");

            } catch (SecurityException se) {
                System.out.println(ColorStatic.RED + "directory " + tempDirect + "could not be created." + ColorStatic.RESET);
            }
        }

        this.structure = struct;
        File file = new File(tempDirect.getAbsoluteFile() + "/temp_" + Long.toString(randID) + ".txt");

        SmarterJmolAdapter adapter = new SmarterJmolAdapter();
        JPanel dummyPanel = new JPanel();
        JmolViewer viewer = JmolViewer.allocateViewer(dummyPanel, adapter);

        try {
            String pdb = structure.toPDB();
            viewer.openStringInline(pdb);
        } catch (Exception e) {
            System.out.println(e.getMessage());

        }

        String command = "select all; isoSurface resolution 0.5"
                + " ignore (solvent) " + typeOfSurface + " "
                + sizeOfSphere.toString() + ";"
                + " y = getProperty(\"isosurfaceData\");"
                + " x = y[\"vertices\"];"
                + " write VAR x " + "\"" + file.getAbsolutePath().replace("\\", "\\\\") + "\";";

        viewer.evalString(command);

        System.out.println(command);
        List<Double[]> vertices = new ArrayList<>();

        for (int i = 200; i > 0; i--) {

            if (i == 1) {

                System.out.println(ColorStatic.GREEN
                        + "Found file " + file.getAbsolutePath() + " - " + file.length() + ColorStatic.RESET);

                FileArrayProvider fap = new FileArrayProvider();

                try {
                    vertices = fap.readLines(file.getAbsolutePath());

                    System.out.println(ColorStatic.GREEN + "vertice dimension "
                            + vertices.size()
                            + ColorStatic.RESET);

                } catch (IOException ex) {
                    System.out.println(ex.getMessage());
                }
                file.deleteOnExit();
                break;

            } else {
                try {
                    Thread.sleep(100);
                } catch (InterruptedException ex) {
                    //do nothing
                }
            }
        }

        viewer.evalString("initialize;");
        viewer.evalString("zap");
        return vertices;
    }

    private boolean isCompletelyWritten(File file) {
        boolean resp = false;
        RandomAccessFile stream = null;

        try {
            stream = new RandomAccessFile(file, "rw");
            resp = true;
            System.out.print("+");
        } catch (Exception e) {
            System.out.print("-");
        } finally {
            if (stream != null) {
                try {
                    stream.close();
                } catch (IOException e) {
                    System.out.println(e.getMessage());
                }
            }
        }
        return resp;
    }

    public Double getResolution() {
        return resolution;
    }

    public void setResolution(Double resolution) {
        this.resolution = resolution;
    }

    public Double getSizeOfSphere() {
        return sizeOfSphere;
    }

    public void setSizeOfSphere(Double sizeOfSphere) {
        this.sizeOfSphere = sizeOfSphere;
    }

    public String getTypeOfSurface() {
        return typeOfSurface;
    }

    public void setTypeOfSurface(String typeOfSurface) {
        this.typeOfSurface = typeOfSurface;
    }

    @Override
    public void run() {

    }

    public Structure getStructure() {
        return structure;
    }

    public void setStructure(Structure structure) {
        this.structure = structure;
    }

    private class MyJmolViewer {

        Logger logger;
        JmolViewer viewer;
        JmolAdapter adapter;
        JmolStatusListener statusListener;

        Structure structure;

        MyJmolViewer() {

            this.logger = LoggerFactory.getLogger(JmolPanel.class);
            this.adapter = new SmarterJmolAdapter();
            viewer = JmolViewer.allocateViewer(this, this.adapter, null, null, null, null, this.statusListener);
        }

        public JmolViewer getViewer() {
            return this.viewer;
        }

        public void executeCmd(String rasmolScript) {
            this.viewer.evalString(rasmolScript);
        }

        public void setStructure(Structure s) {
            if (s == null) {
                return;
            }

            this.structure = s;

            try {
                String pdb = s.toPDB();
                this.viewer.openStringInline(pdb);
            } catch (Exception e) {
                logger.error(e.getMessage());
            }
        }

        public void openStringInline(String pdbFile) {
            this.viewer.openStringInline(pdbFile);
        }
    }

}
