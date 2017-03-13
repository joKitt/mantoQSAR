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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.biojava.nbio.structure.Atom;
import org.mantoQSAR.core.DescriptorSet;
import org.mantoQSAR.core.Molecule;
import org.mantoQSAR.core.util.MoleculeTools;
import org.mantoQSAR.core.util.Projection;


public class ShapeDescriptorGroup extends DescriptorGroup{

    List<Double[]> surfacePoints;
    List<Double> atomMass;
    
    public ShapeDescriptorGroup(){
        super();
    }
    
    
    public ShapeDescriptorGroup(int observationID, int descriptorID) {
        super(observationID, descriptorID);
    
        // number of descriptors in this group
        int nD = 11;

        this.descriptorList = new ArrayList<>();
        for (int i = 0; i < nD; i++) {
            this.descriptorList.add(new Descriptor("placeholder", 0.0));
        };
    }

    @Override
    public void calcDescriptor(){

        Molecule m = this.getMolecule();
        
        surfacePoints = m.getSurface();
        atomMass = MoleculeTools.getAtomMass(m.getAtomList());
        this.descriptorList = this.calcDescriptorValue();
        
} 
    
    private List<Descriptor> calcDescriptorValue(){
        
       Molecule m = this.getMolecule();
       DescriptorSet descriptorSet = this.getDescriptorSet();
       
       List<Descriptor> descrList = new ArrayList<>();

       String string = descriptorSet.getDescriptor().name;
        string = string.substring(0,1).toUpperCase() + string.substring(1);
    
        if(!string.substring(0,1).equalsIgnoreCase("_")){
        string = "_" + string;
        }

        Double[] centerPmass = MoleculeTools.getCenter(m.getAtomList(), atomMass);
        
        List<Double[]> cPList = new ArrayList<>();
        cPList.add(centerPmass);
        Projection p = new Projection();
        List<Double> radi = Arrays.asList(p.getAbsDistance(cPList, surfacePoints));
        
        List<Double> sortRadi = new ArrayList<>(); 
        for (Double dd:radi){
            sortRadi.add(dd); 
        }
        
        Collections.sort(sortRadi); 
        int center = sortRadi.size()/2; 
        double median = sortRadi.get(center);
        double max = Collections.max(radi);
        double min = Collections.min(radi);
        
        double res = Math.pow(descriptorSet.getSurface().getResolution(), 2); 
        
        int nAtom = m.getAtomList().size();
        int nAminoAcid = m.getAsList().size();
        
        double charge = 0.0; 
        double mass = 0.0; 
        
        for(Atom a:m.getAtomList()){
            charge = charge + a.getOccupancy();
            mass = mass + a.getElement().getAtomicMass();
        }
 
        descrList.add(new Descriptor("shapeMax" + string, max/median));
        descrList.add(new Descriptor("shapeMin" + string, min/median));
        
        descrList.add(new Descriptor("shapeFactor" + string, (max-min)/median));
        
        descrList.add(new Descriptor("nSurfP" + string, (double) radi.size()));
        descrList.add(new Descriptor("surfArea" + string, (double) radi.size()*res));

        descrList.add(new Descriptor("nAtom" + string, (double) nAtom));
        descrList.add(new Descriptor("nAAcid" + string, (double) nAminoAcid));
        descrList.add(new Descriptor("mass" + string, mass));

        descrList.add(new Descriptor("dens" + string, mass/(radi.size()*res)));

        descrList.add(new Descriptor("charge" + string, charge));
        descrList.add(new Descriptor("chargeDens" + string, charge/nAtom));

       return descrList; 
    }

    public List<Double[]> getSurfacePoints() {
        return surfacePoints;
    }

    public void setSurfacePoints(List<Double[]> surfacePoints) {
        this.surfacePoints = surfacePoints;
    }

    public List<Double> getAtomMass() {
        return atomMass;
    }

    public void setAtomMass(List<Double> atomMass) {
        this.atomMass = atomMass;
    }
}
