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


package org.mantoQSAR.core;

import java.util.ArrayList;
import java.util.List;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.Structure;
import org.mantoQSAR.core.math.Plane;
import org.mantoQSAR.core.util.Isosurface;
import org.mantoQSAR.core.util.MoleculeTools;
import org.mantoQSAR.core.util.Projection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class Molecule{

Logger logger; 

private Structure structure; 

private Double ionicStrength = 0.0;
private Double dielectricConstant = 2.0;
private List<Double[]> surface;
private final Double resolution = 0.5; 

    public Molecule(){
        
        logger = LoggerFactory.getLogger(Molecule.class);
        this.surface = new ArrayList<>();
    }

    public Structure getStructure() {
        return structure;
    }

    public void setStructure(Structure structure) {
        this.surface.clear();
        this.structure = structure;

    }
    
    public List<Double[]> getResiduePosition(String type) {

        List<Double[]> resPos = new ArrayList<>();
        List<Group> asList = this.getAsList();
        
        switch (type) {
            case "mass":

                for (Group a : asList) {

                    List<Atom> aList = a.getAtoms();
                    resPos.add(MoleculeTools.getCenter(aList, MoleculeTools.getAtomMass(aList)));
                }
                break;
            case "geo":
                for (Group a : asList) {
                    List<Atom> aList = a.getAtoms();
                    resPos.add(MoleculeTools.getCenter(aList));
                }
                break;
            default:
                for (Group a : asList) {
                    List<Atom> aList = a.getAtoms();
                    resPos.add(MoleculeTools.getCenter(aList));
                }
                break;
        }

        return resPos;
    }

    public List<Double> getHydrophobicityConstant() {

        return MoleculeTools.getHydrophobicityConstant(this.getAsList());

    }
    /*
    public void centerMolecule(){

        logger.info("center molecule ");
        
        if(this.structure == null){
            return; 
        }
        
        Double[] center = MoleculeTools.getCenter(this.getAtomList());

        List<Chain> model = this.structure.getModel(0);
       // List<Atom> atomList = new ArrayList<>();
        
        for(int i = 0; i < model.size(); i++){
            Chain c = model.get(i);
            
            List<Group> aG = c.getAtomGroups();
            
            for(int j = 0; j < aG.size(); j++){
                Group g = aG.get(j);
                
                List<Atom> aList = g.getAtoms();
                
                for(int k = 0; k< aList.size(); k++){
                    
                    Atom a = aList.get(k);
                    
                    double[] coord = new double[3];
                    coord[0] = a.getCoords()[0] - center[0];
                    coord[1] = a.getCoords()[1] - center[1];
                    coord[2] = a.getCoords()[2] - center[2];
                    
                    a.setCoords(coord); 
                    aList.set(k, a);
                }
                g.setAtoms(aList);
                aG.set(j, g);
            }
            c.setAtomGroups(aG);
            model.set(i, c);
        }
        
        this.structure.setModel(0, model);
       
    }
    */

     public List<Double> getOccupancy() {
         // pqr files code atom charge as occupancy
         
        List<Double> atomValue = new ArrayList<>();
        List<Atom> aList = this.getAtomList(); 
        
        for (Atom a : aList) {
            atomValue.add(a.getOccupancy());
        }
        return atomValue;
    }
    
    public List<Double[]> getSurface() {

       if(this.surface.isEmpty()){
        this.surface = this.calcSurface();
        }
        return this.surface; 
       
    }
     
    public void clearSurface(){
        
        this.surface.clear();
    }
    
    private List<Double[]> calcSurface(){
        Isosurface isosurface = new Isosurface();
        isosurface.setResolution(this.resolution);

        return isosurface.getSurface(this.getStructure());
    }
    

    public List<Group> getAsList() {
        if(this.structure == null){
            return null; 
        }
        
        List<Chain> model = this.structure.getModel(0);
        List<Group> asList = new ArrayList<>();
        
        for (Chain c : model) {
            List<Group> atomGroup = c.getAtomGroups();

            atomGroup.stream().forEach((a) -> {
                asList.add(a);
            });
        }
        return asList; 
    }

    public List<Double> getAtomMass() {
        return (MoleculeTools.getAtomMass(this.getAtomList()));
    }


    public List<Atom> getAtomList(){
        
        if(this.structure == null){
            return null; 
        }
        
        List<Chain> model = this.structure.getModel(0);
        List<Atom> atomList = new ArrayList<>();
        
        for (Chain c : model) {
            List<Group> atomGroup = c.getAtomGroups();

            for (Group a : atomGroup) {
                List<Atom> aList = a.getAtoms();
                for (Atom atom : aList) {
                    atomList.add(atom);
                }
            }
        }
        return atomList; 
        
    }
    
    public List<Double[]> getAtomPosition() {
        return MoleculeTools.getAtomPosition(this.getAtomList());
    }

    public Double getIonicStrength() {
        return ionicStrength;
    }

    public void setIonicStrength(Double ionicStrength) {
        this.ionicStrength = ionicStrength;
    }

    public Double getDielectricConstant() {
        return dielectricConstant;
    }

    public void setDielectricConstant(Double dielectricConstant) {
        this.dielectricConstant = dielectricConstant;
    }

    public List<Atom> getAtomPatch(Double[] vector, Double size){
          throw new UnsupportedOperationException("Not implemented yet");
    }
    
    public List<Double[]> getSurfacePatch(Double[] vector, Double size) {

        List<Double[]> refP = this.getSurface();
        List<Double[]> surfPart = new ArrayList<>();
        
        // calculate a plane, properties not critical, as plane only needed for reference
        Plane plane = new Plane(this.getAtomList(), vector, 100.0, 5.0, 0.0); 
        List<Double[]> planeP = plane.getPlane();
        Projection p = new Projection(); 
        
        Double[] dist2plane = p.getAbsDistance(planeP, refP);
        p = null; 
       
        Double minDistVal = 10000.0; 
        for (int i = 0; i < dist2plane.length; i++) {
            if(dist2plane[i] < minDistVal){
                minDistVal = dist2plane[i]; 
            }
            if (dist2plane[i] < size) {
                surfPart.add(refP.get(i));
            }
        }
        
        return surfPart;
    }
    
    public boolean[] getSurfacePatchIO(Double[] vector, Double size) {

        List<Double[]> refP = this.getSurface();
        boolean[] surfPartIO = new boolean[refP.size()];
        
        // calculate a plane, properties not critical, as plane only needed for reference
        Plane plane = new Plane(this.getAtomList(), vector, 100.0, 5.0, 0.0); 
        List<Double[]> planeP = plane.getPlane();
       
        Projection p = new Projection(); 
        
        Double[] dist2plane = p.getAbsDistance(planeP, refP);
        p = null; 
       
        for (int i = 0; i < dist2plane.length; i++) {

            if (dist2plane[i] < size) {
                surfPartIO[i] = true; 
            }else{
                surfPartIO[i] = false; 
            }
        }

        return surfPartIO;
    }

}
