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

import java.awt.Color;
import java.awt.Dimension;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.Collections;
import java.util.List;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import org.biojava.nbio.structure.Structure;
import org.jmol.api.JmolViewer;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.Screen;
import org.mantoQSAR.core.descriptor.*;
import org.mantoQSAR.core.math.Matrix;
import org.mantoQSAR.core.math.Vector;
import org.mantoQSAR.core.util.ColorStatic;
import org.mantoQSAR.core.util.MoleculeTools;
import org.mantoQSAR.gui.resources.ColorRange;


public class MoleculeViewPanel extends JPanel implements PropertyChangeListener {

   
    private Screen screen;
    private MoleculeController moleculeControl;
    private String panelTitle;
    private JmolPanel viewPanel; 
    private int styleID = 0;
    private ColorRange colorRange; 
    
    private static final int GAP = 5;

    private boolean showDescriptor = false;

    
    public MoleculeViewPanel() {
        
        moleculeControl = MoleculeController.getInstance();
        moleculeControl.addPropertyChangeListener(this);
        
        screen = moleculeControl.getScreen();
        screen.addPropertyChangeListener(this);

        this.panelTitle = "Molecule";
        this.colorRange = new ColorRange();

        initPanel();

        this.setBackground(Color.WHITE);
        this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
        this.setBorder(BorderFactory.createEmptyBorder(GAP, GAP, GAP, GAP));
        this.add(viewPanel);
    }

    protected final void initListeners() {
        // get screen
        moleculeControl.addPropertyChangeListener(this);
        screen.addPropertyChangeListener(this);
    }

    private void initPanel() {
        
        this.viewPanel = moleculeControl.getViewPanel();
        viewPanel.setBackground(Color.WHITE);
        viewPanel.setBorder(new TitledBorder(this.panelTitle));

        viewPanel.setPreferredSize(new Dimension(500, 500));
        viewPanel.executeCmd("background white; frank off; set antialiasdisplay on");
    }

    public void setShowDescriptor(boolean b) {

        this.showDescriptor = b;

        if (b == true) {
            this.displayDescriptor(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor());
        }
    }

    public void displayDescriptor(int molPos, int descPos) {

        Integer moleculePos = this.moleculeControl.getActivePosition();
        Integer descrPos = this.moleculeControl.getActiveDescriptor();

        executeCmd("color chains;");

        if (screen.getDescriptorGroup(moleculePos, descrPos) instanceof SurfaceDescriptorGroup) {
            System.out.println(ColorStatic.BLUE + "Display surface for surface descriptor" + ColorStatic.RESET);

            executeCmd("draw DELETE;");
            
            SurfaceDescriptorGroup sdG = (SurfaceDescriptorGroup) screen.getDescriptorGroup(moleculePos, descrPos);
            List<Double[]> surf = sdG.getMolecule().getSurface();
            List<Double> valueList = sdG.getSurfaceValue();
            this.displayPoints(surf, null, 0.6);

        }

        // show plane descriptor
        if (screen.getDescriptorGroup(moleculePos, descrPos) instanceof PlaneDescriptorGroup) {
            System.out.println(ColorStatic.BLUE + "Display surface for plane descriptor" + ColorStatic.RESET);

            PlaneDescriptorGroup sdG = (PlaneDescriptorGroup) screen.getDescriptorGroup(moleculePos, descrPos);
            
            if(moleculeControl.getOrientationID() == null){
              moleculeControl.setOrientationID(sdG.getPreferredOrientationIndex());
            }

            List<Double[]> plane = sdG.getPlane(this.moleculeControl.getOrientationID());
            Double[] vector = sdG.getVector().get(this.moleculeControl.getOrientationID());

            List<Double> valueList = sdG.getPlaneValue(this.moleculeControl.getOrientationID());
           
            executeCmd("draw DELETE;");
            executeCmd("set echo top left;");
            executeCmd("color echo gray;");
            executeCmd("font echo 12 arial;");
            executeCmd("echo \"plane " + (this.moleculeControl.getOrientationID()+1) + ": " + String.format("%.2f", vector[0]) 
                                 + " " + String.format("%.2f", vector[1]) 
                                 + " " + String.format("%.2f", vector[2]) 
                            + " \";");

            this.displayPoints(plane, null, 0.6);

        }

        // show patch descriptor
        if (screen.getDescriptorGroup(moleculePos, descrPos) instanceof PatchDescriptorGroup) {
         System.out.println(ColorStatic.BLUE + "Display surface for patch descriptor" + ColorStatic.RESET);

            PatchDescriptorGroup sdG = (PatchDescriptorGroup) screen.getDescriptorGroup(moleculePos, descrPos);
            
            if(moleculeControl.getOrientationID() == null){
              moleculeControl.setOrientationID(sdG.getPreferredOrientationIndex());
            }

            executeCmd("draw DELETE;");
            List<Double[]> surfPatch = sdG.getPatch(this.moleculeControl.getOrientationID());
            List<Double> valueList = sdG.getPatchValue(this.moleculeControl.getOrientationID());
            Double[] vector = sdG.getVector().get(this.moleculeControl.getOrientationID());
            executeCmd("set echo top left;");
            executeCmd("color echo gray;");
            executeCmd("font echo 12 arial;");
            executeCmd("echo \"surface patch " + (this.moleculeControl.getOrientationID()+1) + ": " + String.format("%.2f", vector[0]) 
                                 + " " + String.format("%.2f", vector[1]) 
                                 + " " + String.format("%.2f", vector[2]) 
                            + " \";");

            this.displayPoints(surfPatch, null, 0.6);
        } 

        // show descriptor for sphere descriptor
        if (screen.getDescriptorGroup(moleculePos, descrPos) instanceof SphereDescriptorGroup) {
            System.out.println(ColorStatic.PURPLE + "Display surface for sphere descriptor" + ColorStatic.RESET);

            SphereDescriptorGroup sdG = (SphereDescriptorGroup) screen.getDescriptorGroup(moleculePos, descrPos);
            
            // delete any previous descriptors displayed
            executeCmd("draw DELETE;");
            
            List<Double[]> surf = sdG.getMolecule().getSurface();
            this.displayPoints(surf, null, 0.6);

        }

        if (screen.getDescriptorGroup(moleculePos, descrPos) instanceof ShapeDescriptorGroup) {
            System.out.println(ColorStatic.PURPLE + "Display surface for shape descriptor" + ColorStatic.RESET);

            ShapeDescriptorGroup sdG = (ShapeDescriptorGroup) screen.getDescriptorGroup(moleculePos, descrPos);
            
            // delete any previous descriptors displayed
            executeCmd("draw DELETE;");
            
            List<Double[]> surf = sdG.getMolecule().getSurface();

            this.displayPoints(surf, null, 0.6);

        } 
        
        String cmdString = "write POVRAY descr_" + Integer.toString(molPos) +  "_" 
                            + Integer.toString(descrPos) + "_"
                            + moleculeControl.getOrientationID();
        executeCmd(cmdString);
        

    }

    public void displayVector(List<Double[]> pointA, Double[] pointB) {

        for (int i = 0; i < pointA.size(); i++) {
            Double[] pA = pointA.get(i);

            viewPanel.executeCmd("draw line" + Integer.toString(i) + " {"
                    + pA[0].toString() + " "
                    + pA[1].toString() + " "
                    + pA[2].toString()
                    + "} {"
                    + pointB[0].toString() + " "
                    + pointB[1].toString() + " "
                    + pointB[2].toString()
                    + "}; color $line" + Integer.toString(i) + " red;");
        }
    }

    public void displayPoints(List<Double[]> pointMatrix, List<Double> value, Double size) {

        double pointDiameter = 0.4;
        if(size != null){
            pointDiameter = size;
        }

        if(value!= null){
            this.colorRange = new ColorRange();
            this.colorRange.setSigmoidal(false);
            this.colorRange.setScheme(2);

     
            if(Collections.max(value) > 0){
            this.colorRange.setMaxPos(Collections.max(value));
           
            }
            
            if(Collections.min(value) < 0){
             this.colorRange.setMin(Collections.min(value));
       
            }
        }

        int n = pointMatrix.size();

        if (value != null && value.size() == n) {
            for (int i = 0; i < n; i++) {
             Color color = this.colorRange.getColor(value.get(i));
                
                Double[] coord = pointMatrix.get(i);
                String cmdString = "draw c" + Integer.toString(i) + " circle {"
                        + coord[0].toString() + " "
                        + coord[1].toString() + " "
                        + coord[2].toString() + "} diameter " + Double.toString(pointDiameter)
                        + " color {" + color.getRed() + " " + color.getGreen() + " " + color.getBlue() + "};";
                
                viewPanel.executeCmd(cmdString);
             }
        } else {

            for (int i = 0; i < n; i++) {
                Double[] coord = pointMatrix.get(i);
                String cmdString = "draw c" + Integer.toString(i) + " circle {"
                        + coord[0].toString() + " "
                        + coord[1].toString() + " "
                        + coord[2].toString() + "} diameter " + Double.toString(pointDiameter)
                          + " color {0.84,0.176, 0.125};";
                viewPanel.executeCmd(cmdString);
            }
        }
    }

    public void displayMolecule(int moleculePos) {

        String cmdString; 
        // display molecule in screen position i
        Molecule m = null;
        try {
            m = screen.getMolecule(moleculePos);
            System.out.println(ColorStatic.PURPLE + "show molecule with " + m.getAtomList().size() + " atoms"
                    + ColorStatic.RESET);
        
        } catch (NullPointerException e) {
            System.out.println(ColorStatic.PURPLE + e.getMessage() + ColorStatic.RESET);
        }
        
        if (m != null) {
            this.panelTitle = screen.getObservationSetList().get(moleculePos).getName();
            this.displayStructure(m.getStructure());

            // draw semi transperent SASA surface
            executeCmd("select all; isosurface surf1 ignore(solvent) sasurface;");
            executeCmd("color isosurface white translucent 0.50;");
            
            Double[] center = MoleculeTools.getCenter(m.getAtomList());

            if(center != null){
                cmdString = "centerAt absolute {" + center[0].toString() + " " + center[1].toString() + "" + center[2].toString() + "};";
                executeCmd(cmdString);
            }
        }
        
    }

    public void displayStructure(Structure structure) {

        String pdb = structure.toPDB();
        executeCmd("zap;");
        this.viewPanel.getViewer().openStringInline(pdb);
        this.viewPanel.getViewer().evalString("select *;");
        this.setStyle(this.styleID);

    }

    public void displayStructureFile(File f) {


        executeCmd("zap;");
        System.out.println("Open file: " + f.toString());
        viewPanel.getViewer().openFile(f.toString());
        setStyle(this.styleID);
    }

    public void executeCmd(String rasmolScript) {
        viewPanel.getViewer().evalString(rasmolScript);
    }

    public JmolPanel getJmolPanel() {
        return viewPanel;
    }

    public JmolViewer getJmolViewer() {
        return viewPanel.getViewer();
    }

    
    public ColorRange getColorRange() {
        return colorRange;
    }

    public void setColorRange(ColorRange colorRange) {
        this.colorRange = colorRange;
    }

    private void setStyle(int styleID) {

        switch (styleID) {
            
// cartoon style
            case 0:
               executeCmd("select all;"         
                        + "background WHITE; "
                        + "frank off; "
                        + "hbonds off; "
                        + "wireframe off; "
                        + "spacefill off; "
                        + "trace off; "
                        + "set ambient 40; "
                        + "set specpower 40; "
                        + "cartoon only; "
                        + "color group; "
                        + "spin off; ");
                break;

            // gray wireframe
            case 1:

               executeCmd("select all;"   
                        + "background WHITE; "
                        + "frank off; "
                        + "hbonds off; "
                        + "wireframe off; "
                        + "spacefill off; "
                        + "trace off; "
                        + "set ambient 40; "
                        + "set specpower 40; "
                       + " set diffuse 100; "
                        + "wireframe 0.3; "
                        + "color [224, 224, 224]; ");
                break;

            // default is comic style
            default:
                executeCmd("select all;"
                        + "background WHITE; "
                        + "frank off; "
                        + "hbonds off; "
                        + "wireframe off; "
                        + "spacefill off; "
                        + "trace off; "
                        + "set ambient 40; "
                        + "set specpower 40; "
                        + "cartoon only; "
                        + "color group; "
                        + "spin off; ");
        }
    }

    @Override
    public void propertyChange(PropertyChangeEvent evt) {

        if (evt.getPropertyName().equalsIgnoreCase(MoleculeGuiStatic.CHANGE_ACTIVE_OBSERVATION)) {

            this.displayMolecule(this.moleculeControl.getActivePosition());
            if (this.showDescriptor == true) {
                
                // reset orientation ID to get prefered orientation for new molecule to display
                this.moleculeControl.setOrientationID(null);
                this.displayDescriptor(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor());
            }
        }

        if (evt.getPropertyName().equalsIgnoreCase(MoleculeStatic.CHANGE_LOAD_PROPERTY_SETTING)) {
            this.displayMolecule(this.moleculeControl.getActivePosition());
        }

        
        if (evt.getPropertyName().equalsIgnoreCase(MoleculeGuiStatic.CHANGE_ACTIVE_DESCRIPTOR)) {
            System.out.println(ColorStatic.PURPLE + " fired change_active_descriptor to new descriptor "
                    + this.moleculeControl.getActiveDescriptor().toString() + ColorStatic.RESET);

            this.setShowDescriptor(true);
            this.displayDescriptor(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor());
        }
        
         if (evt.getPropertyName().equalsIgnoreCase(MoleculeGuiStatic.CHANGE_ACTIVE_ORIENTATION)) {
            System.out.println(ColorStatic.PURPLE + " fired change_active_descriptor to new descriptor "
                    + this.moleculeControl.getActiveDescriptor().toString() + ColorStatic.RESET);

            this.setShowDescriptor(true);
            this.displayDescriptor(this.moleculeControl.getActivePosition(), this.moleculeControl.getActiveDescriptor());
        }
    }

    public int getStyleID() {
        return styleID;
    }

    public void setStyleID(int styleID) {
        this.styleID = styleID;
        this.displayMolecule(this.moleculeControl.getActivePosition());
    }

    /*  Rotates the molecular view, for surface to be at the bottom. 
        Rotation matrices in Jmol are indicated either as a 3-row, 
       3-column array or by 9 numbers set off by double brackets.
       */
    public void rotate(Double[] vector){
        
       System.out.println(ColorStatic.PURPLE + "rotate called on structure" + ColorStatic.RESET);
    
        Double[] vectorTarget = {0.0, -1.0,  0.0}; 

        // norm vector before use
        for (int i = 0; i < vector.length; i++) {
            vector[i] = vector[i] / Vector.norm(vector);
            vectorTarget[i] = vectorTarget[i] / Vector.norm(vectorTarget);
        }

        // calculate cross product of vectors
        Double[] n = Vector.crossProduct(vector, vectorTarget);
        for (int i = 0; i < n.length; i++) {
            n[i] = n[i] / Vector.norm(n);
        }
        Double cosAlpha = vector[0]*vectorTarget[0]+ vector[1]*vectorTarget[1]+ vector[2]*vectorTarget[2];
        Double sinAlpha = Math.sqrt(Math.pow(n[0], 2) + Math.pow(n[1], 2) +Math.pow(n[2], 2)) ;
        
        Matrix RAlpha = new Matrix(3,3);
        RAlpha.data[0][0] = Math.pow(n[0],2) * (1-cosAlpha) + cosAlpha;
        RAlpha.data[1][0] = (n[1]*n[0]*(1-cosAlpha)) + (n[2]*sinAlpha);
        RAlpha.data[2][0] = (n[2]*n[0]*(1-cosAlpha)) - (n[1]*sinAlpha);
        
        RAlpha.data[0][1] = (n[0]*n[1]*(1-cosAlpha)) - (n[2]*sinAlpha);
        RAlpha.data[1][1] = Math.pow(n[1],2) * (1-cosAlpha) + cosAlpha;
        RAlpha.data[2][1] = (n[2]*n[1]*(1-cosAlpha)) + (n[0]*sinAlpha);
        
        
        RAlpha.data[0][2] = (n[0]*n[2]*(1-cosAlpha)) + (n[1]*sinAlpha);
        RAlpha.data[1][2] = (n[1]*n[2]*(1-cosAlpha)) - (n[0]*sinAlpha);
        RAlpha.data[2][2] = Math.pow(n[2],2) * (1-cosAlpha) + cosAlpha;
        
        Double angleRad = Math.acos(cosAlpha); // rad 
        Double angleDegree = angleRad * (180/Math.PI);
        System.out.println(ColorStatic.GREEN + "angle to rotate by is " + angleDegree + " degree" + ColorStatic.RESET);

        String cmdString = "rm = [["   + RAlpha.data[0][0] + " "
                                       + RAlpha.data[0][1] + " "
                                       + RAlpha.data[0][2] + " "
                                       + RAlpha.data[1][0] + " "
                                       + RAlpha.data[1][1] + " "
                                       + RAlpha.data[1][2] + " "
                                       + RAlpha.data[2][0] + " "
                                       + RAlpha.data[2][1] + " "
                                       + RAlpha.data[2][2] + " "
                                       + "]];";        
        executeCmd(cmdString);
        executeCmd("rotate @rm");
        executeCmd("select all; center selected;");
    }
    
    public void colorAminoAcids(int index){
        
        switch(index){
        case 1:
            executeCmd("select[lys]; color red;");
            executeCmd("select[arg]; color green;");
            break; 
        case 2:
            executeCmd("select[lys]; color {150 150 150};");
            executeCmd("select[arg]; color {180 180 180};");
            break;
        case 3:
            executeCmd("select[lys]; color red;");
            executeCmd("select[arg]; color red;");
            break;     
        default: 
    }
    }

    private void labelAminoAcids() {
        executeCmd("select[lys]*.CA; label %n %r; "
                    + " label %n %r;"
                    + " font label 12 arial; color label black;"
                    + " set labelalignment right; set labeloffset 20 -20; "
                    + " set labelfront on;");
        
        executeCmd("select[arg]*.CA; "
                    + " label %n %r;"
                    + " font label 12 arial; color label black;"
                    + " set labelalignment right; set labeloffset 20 -20; "
                    + " set labelfront on;");
    }
} 
   
