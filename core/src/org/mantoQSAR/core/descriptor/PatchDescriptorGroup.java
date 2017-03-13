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
import org.mantoQSAR.core.DescriptorSet;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.util.ColorStatic;
import org.mantoQSAR.core.util.Projection;


public class PatchDescriptorGroup extends OrientationDescriptorGroup {

    List<Double> valueMap;

    public PatchDescriptorGroup(){
        super();
    }
    
    
    public PatchDescriptorGroup(int observationID, int descriptorID) {
        super(observationID, descriptorID);

        this.valueMap = new ArrayList<>();

        int nD = 28;

        this.descriptorList = new ArrayList<>();
        for (int i = 0; i < nD; i++) {
            this.descriptorList.add(new Descriptor("placeholder", 0.0));
        };


    }

    @Override
    public void calcDescriptor() {

        Molecule m = this.getMolecule();
        DescriptorSet descriptorSet = this.getDescriptorSet();        
        List<Double[]> surfacePoints = m.getSurface();

        this.valueMap = new ArrayList<>(surfacePoints.size());

        Projection p = new Projection();
        
        switch (descriptorSet.getSurface().getProperty()) {
            case "esp":
                valueMap = p.calcProjection(surfacePoints,
                        m.getAtomPosition(),
                        m.getOccupancy(),
                        descriptorSet.getSurface().getMapFunction(),
                        0.0);
                break;

            case "hyd":
                valueMap = p.calcProjection(surfacePoints,
                        m.getResiduePosition("mass"),
                        m.getHydrophobicityConstant(),
                        descriptorSet.getSurface().getMapFunction(),
                        0.0);

                break;

            default:
               
                System.out.println(ColorStatic.RED + "Property to be mapped in descriptor calculation not identified." + ColorStatic.RESET);
                return;
        }
        
        p = null; 

        for (Double[] vector1 : this.vector) {

            boolean[] patchPointIO = m.getSurfacePatchIO(vector1, descriptorSet.getProjection().getSize());

            int patchCount = 0; 
            for(int i = 0; i < patchPointIO.length; i++){
                if(patchPointIO[i] == true){
                    patchCount ++; 
                }
            }

            List<Double> valueMap2 = new ArrayList<>(patchCount); 
            for(int i = 0; i < patchPointIO.length; i++){
                if(patchPointIO[i] == true){
                    valueMap2.add(valueMap.get(i));
                }   
            }
 
            this.descriptorDetail.add((ArrayList<Descriptor>) this.calcDescriptorValue(valueMap2));

        }

        logger.info("Select orientation based on " + descriptorSet.getProjection().getSelectID()
                + " selected for " + descriptorSet.getProjection().getSelectFunction());

        int descrInd = this.getDescriptorPosition(descriptorSet.getProjection().getSelectID());
        logger.info("Descriptor position in set is " + descrInd);

        int prefInd = this.getOrientation(descrInd, descriptorSet.getProjection().getSelectFunction());
        logger.info("Prefered orientation identified with " + descrInd);

        this.descriptorList = this.descriptorDetail.get(prefInd);

        this.setCALC_STATE(true);
        
    }

    private List<Descriptor> calcDescriptorValue(List<Double> valueMap) {

        DescriptorSet descriptorSet = this.getDescriptorSet();
        List<Descriptor> descrList = new ArrayList<>();

        String string = descriptorSet.getDescriptor().name;
        string = string.substring(0, 1).toUpperCase() + string.substring(1);

        if (!string.substring(0, 1).equalsIgnoreCase("_")) {
            string = "_" + string;
        }

        double size = (double) valueMap.size();
        double res = Math.pow(descriptorSet.getSurface().getResolution(), 2);

        double c = size / res;
        descrList.add(new Descriptor("totalSurf" + string, c));
        descrList.add(new Descriptor("nSurfP" + string, size));

        double sum = 0.0;
        for (Double dd : valueMap) {
            sum = sum + dd;
        }

        double mean = sum / size;

        descrList.add(new Descriptor(("sumSurf" + string), sum));
        descrList.add(new Descriptor(("mean" + string), mean));

        List<Double> sortValueMap = new ArrayList<>();
        for (Double dd : valueMap) {
            sortValueMap.add(dd);
        }

        Collections.sort(sortValueMap);

        int center = sortValueMap.size() / 2;

        double median = sortValueMap.get(center);

        double max = Collections.max(valueMap);
        double min = Collections.min(valueMap);

        descrList.add(new Descriptor(("medianSurf" + string), median));
        descrList.add(new Descriptor(("meanRes" + string), (mean) / res));

        descrList.add(new Descriptor(("max" + string), max));
        descrList.add(new Descriptor(("min" + string), min));

        descrList.add(new Descriptor(("devA" + string), (max - min) / median));
        descrList.add(new Descriptor(("devB" + string), Collections.max(valueMap)));

        List<Double> posValueMap = new ArrayList<>();
        List<Double> negValueMap = new ArrayList<>();
        double sumPosVal = 0.0;
        double sumNegVal = 0.0;

        for (Double dd : valueMap) {
            if (dd > 0.0) {
                posValueMap.add(dd);
                sumPosVal = sumPosVal + dd;
            }
            if (dd < 0.0) {
                negValueMap.add(dd);
                sumNegVal = sumNegVal + dd;
            }
        }

        descrList.add(new Descriptor(("nPos" + string), (double) posValueMap.size()));
        descrList.add(new Descriptor(("nNeg" + string), (double) posValueMap.size()));

        descrList.add(new Descriptor(("relPos" + string), (double) posValueMap.size() / valueMap.size()));
        descrList.add(new Descriptor(("relNeg" + string), (double) posValueMap.size() / valueMap.size()));

        descrList.add(new Descriptor(("sumPos" + string), sumPosVal));
        descrList.add(new Descriptor(("sumNeg" + string), sumNegVal));

        descrList.add(new Descriptor(("averPos" + string), sumPosVal / posValueMap.size()));
        descrList.add(new Descriptor(("averNeg" + string), sumNegVal / negValueMap.size()));

        double binScale = descriptorSet.getDescriptor().binScale;
        int nBin = 10;
        List<Double> binVal = new ArrayList<>();

        for (int i = 1; i <= nBin; i++) {
            double binLow = (-(binScale * nBin) * 0.5) + (i) * binScale;
            double binHigh = (-(binScale * nBin) * 0.5) + (i + 1) * binScale;

            double a = 0.0;
            for (Double valueMap1 : valueMap) {
                if (valueMap1 >= binLow && valueMap1 < binHigh) {
                    a = a + 1;
                }
            }
            binVal.add(a);
        }

        for (int i = 0; i < nBin; i++) {
            descrList.add(new Descriptor(("binAbs" + string + "_" + i), (double) binVal.get(i) / res));
        }

        return descrList;
    }

    @JSON(include=false)
    public List<Double[]> getPatch(int i){
        
        return getPatch(this.vector.get(i));
    }
    
    @JSON(include=false)
    public List<Double[]> getPatch(Double[] vector1){
        Molecule m = this.getMolecule();
        DescriptorSet descriptorSet = this.getDescriptorSet();
        
        return m.getSurfacePatch(vector1, descriptorSet.getProjection().getSize());
    }
    
    @JSON(include=false)
    public List<Double> getPatchValue(int index){
        
        Molecule m = this.getMolecule();
        boolean[] patchPointIO = m.getSurfacePatchIO(this.vector.get(index), this.getDescriptorSet().getProjection().getSize());

            // get number of positives (points in patch)
            int patchCount = 0; 
            for(int i = 0; i < patchPointIO.length; i++){
                if(patchPointIO[i] == true){
                    patchCount ++; 
                }
            }

            List<Double> valueMap2 = new ArrayList<>(patchCount); 
            for(int i = 0; i < patchPointIO.length; i++){
                if(patchPointIO[i] == true){
                    valueMap2.add(this.valueMap.get(i));
                }   
            }
        
        return valueMap2; 
        
    }

    public List<Double> getValueMap() {
        return valueMap;
    }

    public void setValueMap(List<Double> valueMap) {
        this.valueMap = valueMap;
    }  
}
