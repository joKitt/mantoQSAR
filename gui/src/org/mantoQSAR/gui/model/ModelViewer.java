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
package org.mantoQSAR.gui.model;

import java.awt.BorderLayout;
import java.awt.Color;
import javax.swing.BoxLayout;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

public class ModelViewer extends JPanel {

    ModelResultGraph modelGraph;

    public ModelViewer() {

        JPanel content = new JPanel();
        content.setBackground(Color.WHITE);

        content.setLayout(new BoxLayout(content, BoxLayout.Y_AXIS));
        modelGraph = new ModelResultGraph();
        content.add(modelGraph);
        content.add(new ModelResultGraphSetting(modelGraph));

        this.setLayout(new BorderLayout());
        this.add(new JScrollPane(content));
    }
}
