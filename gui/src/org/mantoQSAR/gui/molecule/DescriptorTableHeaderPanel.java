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
import java.awt.FlowLayout;
import javax.swing.JButton;
import javax.swing.JPanel;
import org.mantoQSAR.core.descriptor.MoleculeStatic;


public class DescriptorTableHeaderPanel extends JPanel{
    
    private final JButton btnExportExcel;
    private final MoleculeController moleculeControl; 
    
    public DescriptorTableHeaderPanel() {
        
        moleculeControl = MoleculeController.getInstance(); 
        
        JButton btnUpdateTable = new JButton("Update table");
        btnUpdateTable.addActionListener(moleculeControl);
        btnUpdateTable.setActionCommand(MoleculeStatic.CHANGE_CALC_DESCRIPTOR);
        btnUpdateTable.setContentAreaFilled(false);
        btnUpdateTable.setForeground(Color.LIGHT_GRAY);
                
        btnExportExcel = new JButton("Export2Excel");
        btnExportExcel.addActionListener(moleculeControl);
        btnExportExcel.setActionCommand(MoleculeGuiStatic.SAVE_DESCRIPTOR_TABLE);
        btnExportExcel.setContentAreaFilled(false);
        btnExportExcel.setForeground(Color.BLACK);
        
        this.setLayout(new FlowLayout(FlowLayout.RIGHT));
        this.setBackground(Color.WHITE);
        
        this.add(btnUpdateTable);
        this.add(btnExportExcel);
    }

     
}
