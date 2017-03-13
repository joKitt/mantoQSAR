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

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JPanel;
import org.biojava.nbio.structure.Structure;
import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolStatusListener;
import org.jmol.api.JmolViewer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class JmolPanel extends JPanel implements ActionListener {

    Logger logger;
    JmolViewer viewer;
    JmolAdapter adapter;
    JmolStatusListener statusListener;

    Structure structure;
    private Dimension currentSize;
    private Rectangle rectClip;

    JmolPanel() {

        this.logger = LoggerFactory.getLogger(JmolPanel.class);
        this.adapter = new SmarterJmolAdapter();
        viewer = JmolViewer.allocateViewer(this, this.adapter, null, null, null, null, this.statusListener);

    }

    @Override
    @SuppressWarnings("deprecation")
    public void paint(Graphics g) {
        getSize(this.currentSize);
        g.getClipBounds(this.rectClip);
        this.viewer.renderScreenImage(g, this.currentSize, this.rectClip);
    }

    public JmolViewer getViewer() {
        return this.viewer;
    }

    public JmolAdapter getAdapter() {
        return this.adapter;
    }

    public JmolStatusListener getStatusListener() {
        return this.statusListener;
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

    public void evalString(String rasmolScript) {
        this.viewer.evalString(rasmolScript);
    }

    public void openStringInline(String pdbFile) {
        this.viewer.openStringInline(pdbFile);
    }

    @Override
    public void actionPerformed(ActionEvent ae) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
