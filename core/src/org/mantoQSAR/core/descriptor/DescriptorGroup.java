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
import java.util.List;
import org.biojava.nbio.structure.Structure;

import org.mantoQSAR.core.DescriptorSet;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.ObservationSet;
import org.mantoQSAR.core.Screen;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class DescriptorGroup {

    Logger logger;
    int descriptorSetID = -1; 
    int observationID = -1; 
    
    private boolean CALC_STATE = false; 
    
    public List<Descriptor> descriptorList;
    List<String> name;
    String stIdent;
    List<Double[]> vector = null;

    ArrayList<ArrayList<Descriptor>> descriptorDetail;
    int descriptorNumb;

    public DescriptorGroup(){
       this(0, 0, null); 
    }
    
    
    public DescriptorGroup(int observationID, int descriptorID, List<Double[]> vec){
        this.logger = LoggerFactory.getLogger(DescriptorGroup.class);
 
        this.observationID = observationID;
        this.descriptorSetID = descriptorID; 
        
        this.descriptorList = new ArrayList<>();
        this.descriptorDetail = new ArrayList<>();
        this.stIdent = this.getDescriptorSet().getDescriptor().name;

        if (vec == null) {
            vec = new ArrayList<>();
            Double[] a = {0.0, 0.0, 1.0};
            vec.add(a);
        }
        this.vector = vec;
    }

    public DescriptorGroup(int observationID, int descriptorID){
       this(observationID, descriptorID, null);
    }

    public boolean isCALC_STATE() {
        return CALC_STATE;
    }

    public void setCALC_STATE(boolean CALC_STATE) {
        this.CALC_STATE = CALC_STATE;
    }

    public void setDescriptorDetail(ArrayList<ArrayList<Descriptor>> descriptorDetail) {
        this.descriptorDetail = descriptorDetail;
    }

    
    @JSON(include=false)
    public DescriptorSet getDescriptorSet(){
        return Screen.getInstance().getDescriptorSet(this.descriptorSetID);
    }

    
    @JSON(include=false)
    public Molecule getMolecule(){
         return Screen.getInstance().getMolecule(this.observationID);
    }
   
    @JSON(include=false)
    public ObservationSet getObservationSet(){
        return Screen.getInstance().getObservationSet(this.observationID);
        
     }
    
    public List<Descriptor> getDescriptor() {

        if (this.descriptorList.isEmpty()) {

            this.descriptorDetail.clear();
            this.descriptorList.clear();
            
            this.calcDescriptor();
        }

        return this.descriptorList;
    }

    @JSON(include=false)
    public Structure getStructure() {
        return this.getMolecule().getStructure();
    }



    public List<String> getName() {
        return name;
    }

    public void setName(List<String> name) {
        this.name = name;
    }

    public List<Double[]> getVector() {
        return vector;
    }

    public void setVector(List<Double[]> vector) {
        this.vector = vector;
    }

    public void calcDescriptor() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public ArrayList<ArrayList<Descriptor>> getDescriptorDetail() {
        return descriptorDetail;
    }

    public int getDescriptorSetID() {
        return descriptorSetID;
    }

    public void setDescriptorSetID(int descriptorSetID) {
        this.descriptorSetID = descriptorSetID;
    }

    public int getObservationID() {
        return observationID;
    }

    public void setObservationID(int observationID) {
        this.observationID = observationID;
    }

    public List<Descriptor> getDescriptorList() {
        return descriptorList;
    }

    public void setDescriptorList(List<Descriptor> descriptorList) {
        this.descriptorList = descriptorList;
    }

    public String getStIdent() {
        return stIdent;
    }

    public void setStIdent(String stIdent) {
        this.stIdent = stIdent;
    }

    public int getDescriptorNumb() {
        return descriptorNumb;
    }

    public void setDescriptorNumb(int descriptorNumb) {
        this.descriptorNumb = descriptorNumb;
    }
}
