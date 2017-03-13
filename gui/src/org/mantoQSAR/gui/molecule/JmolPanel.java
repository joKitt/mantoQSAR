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

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import javax.swing.JPanel;
import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolViewer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class JmolPanel extends JPanel{
    private JmolViewer viewer;
    private Logger logger; 
    
    public JmolPanel() {
        
            logger = LoggerFactory.getLogger(JmolPanel.class);
            viewer = JmolViewer.allocateViewer(this, new SmarterJmolAdapter(),
                    null, null, null, null, null);
    }

        public JmolViewer getViewer() {
            return viewer;
        }

        public void executeCmd(String rasmolScript) {
            viewer.evalString(rasmolScript);
        }

        final Dimension currentSize = new Dimension();
        final Rectangle rectClip = new Rectangle();

        @Override
        public void paint(Graphics g) {
            getSize(currentSize);
            Rectangle r = g.getClipBounds(rectClip);
            try {
                viewer.renderScreenImage(g, r.width, r.height);
            } catch (Exception ex) {
                logger.error(ex.getMessage());
            }
        }
    }

