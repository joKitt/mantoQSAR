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
import org.openide.awt.ActionID;
import org.openide.awt.ActionReference;
import org.openide.awt.ActionRegistration;
import org.openide.util.NbBundle.Messages;

@ActionID(
        category = "File",
        id = "org.mantoQSAR.gui.molecule.OpenProject"
)
@ActionRegistration(
        displayName = "#CTL_OpenProject"
)
@ActionReference(path = "Menu/File", position = 1950, separatorBefore = 1875, separatorAfter = 2025)
@Messages("CTL_OpenProject=Open Project")

public final class OpenProject implements ActionListener {

MoleculeController moleculeController = MoleculeController.getInstance();

@Override
    public void actionPerformed(ActionEvent e) {

        moleculeController.openProjectDialog(); 
            }
}
