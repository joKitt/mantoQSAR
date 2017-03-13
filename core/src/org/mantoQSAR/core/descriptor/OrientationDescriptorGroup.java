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


package org.mantoQSAR.core.descriptor;

import flexjson.JSON;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Group;
import org.mantoQSAR.core.math.Vector;
import org.mantoQSAR.core.util.ColorStatic;


public abstract class OrientationDescriptorGroup extends DescriptorGroup{

    List<Atom> atomList;
    List<Group> asList;
    private Integer orientationIndex = null;
    
   public OrientationDescriptorGroup(){
       super();
   }
    
    
    public OrientationDescriptorGroup(int observationID, int descriptorID) {
        super(observationID, descriptorID);
        
        this.atomList = new ArrayList<>();
        this.asList = new ArrayList<>();
        List<Double[]> sphereP = Vector.calcSphere(this.getDescriptorSet().getProjection().orientation);
        this.vector = sphereP; 
    }
    
    @Override
    public List<Descriptor> getDescriptor(){
        
        return descriptorList; 
    }

    public void calcDescriptor(List<Double[]> vec){

        this.vector = vec; 
        this.calcDescriptor();
    }

    public void displayDescriptorValue(int descrNumb){
        List<Double> valueList = new ArrayList<>();
        
        this.descriptorDetail.stream().forEach((dL) -> {
            valueList.add(dL.get(descrNumb).getValue());
        });
        
        for(int i = 0; i < valueList.size(); i++){
            System.out.println(ColorStatic.PURPLE + i + " " + valueList.get(i) + ColorStatic.RESET);
        }
        
    }
    
    @JSON(include=false)
    public int getOrientation(int descrNumb, String selectMode){
        
        List<Double> valueList = new ArrayList<>();
        
        this.descriptorDetail.stream().forEach((dL) -> {
            valueList.add(dL.get(descrNumb).getValue());
        });
       
       double v = 0.0; 
        switch(selectMode){
            case "max":
                 v = Collections.max(valueList);
                break; 
            case "min":
                v = Collections.min(valueList);
                break;
    
            default:
                System.out.println("Invalid orientation selector. Please use max or min.");
                break;
        }
                this.orientationIndex = valueList.indexOf(v);

                return this.orientationIndex;
    }

    @JSON(include=false)
    public int getDescriptorPosition(String selectID) {

       List<Descriptor> dL = this.descriptorDetail.get(0);
       
       List<Integer> pos = new ArrayList<>(); 

       if(selectID == null){
       System.out.println("selectID not set appropriately");
       }
       
       for (int i =0; i < dL.size(); i++){
           Descriptor d = dL.get(i);
       
           if (d != null){
           if(d.name.contains(selectID)){
               pos.add(i); 
           }}
       }
       
       if(pos.isEmpty()){
           System.out.println("No matching descriptor found. To identifier " + selectID);
               }
    
       if(pos.size() > 1){
           System.out.println("Non unique descriptor identifier \"" + selectID + "\" selected. ");
           
           for(int ii:pos){
            System.out.print(dL.get(ii).getName() + ", ");
           }
            System.out.print(" match identifier. Please specify.");
       }
       
    return pos.get(0);
    }

    public List<Atom> getAtomList() {
        return atomList;
    }

    public void setAtomList(List<Atom> atomList) {
        this.atomList = atomList;
    }

    public List<Group> getAsList() {
        return asList;
    }

    public void setAsList(List<Group> asList) {
        this.asList = asList;
    }

    public Integer getOrientationIndex() {
        return orientationIndex;
    }

    public void setOrientationIndex(Integer orientationIndex) {
        this.orientationIndex = orientationIndex;
    }

    public Integer getPreferredOrientationIndex() {
        return orientationIndex;
    }
    
    
}
