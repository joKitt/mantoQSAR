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
package org.mantoQSAR.gui;

import java.awt.Color;
import javax.swing.JTextPane;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;
import javax.swing.text.StyleContext;

public class ColorTextPane extends JTextPane {

    private final int maxChar = 1000000;

    static final Color D_Black = Color.getHSBColor(0.000f, 0.000f, 0.000f);
    static final Color D_Red = Color.getHSBColor(0.000f, 1.000f, 0.502f);
    static final Color D_Blue = Color.getHSBColor(0.667f, 1.000f, 0.502f);
    static final Color D_Magenta = Color.getHSBColor(0.833f, 1.000f, 0.502f);
    static final Color D_Green = Color.getHSBColor(0.333f, 1.000f, 0.502f);
    static final Color D_Yellow = Color.getHSBColor(0.167f, 1.000f, 0.502f);
    static final Color D_Cyan = Color.getHSBColor(0.500f, 1.000f, 0.502f);
    static final Color D_White = Color.getHSBColor(0.000f, 0.000f, 0.753f);
    static final Color B_Black = Color.getHSBColor(0.000f, 0.000f, 0.502f);
    static final Color B_Red = Color.getHSBColor(0.000f, 1.000f, 1.000f);
    static final Color B_Blue = Color.getHSBColor(0.667f, 1.000f, 1.000f);
    static final Color B_Magenta = Color.getHSBColor(0.833f, 1.000f, 1.000f);
    static final Color B_Green = Color.getHSBColor(0.333f, 1.000f, 1.000f);
    static final Color B_Yellow = Color.getHSBColor(0.167f, 1.000f, 1.000f);
    static final Color B_Cyan = Color.getHSBColor(0.500f, 1.000f, 1.000f);
    static final Color B_White = Color.getHSBColor(0.000f, 0.000f, 1.000f);
    static final Color cReset = Color.GRAY;
    static Color colorCurrent = cReset;
    String remaining = "";

    int count = 0;

    public void append(Color c, String s) {
        StyleContext sc = StyleContext.getDefaultStyleContext();
        AttributeSet aset = sc.addAttribute(SimpleAttributeSet.EMPTY, StyleConstants.Foreground, c);
        int len = getDocument().getLength();
        setCaretPosition(len);
        setCharacterAttributes(aset, false);
        replaceSelection(s);

        count = count + s.length();

        if (count > maxChar) {
            try {
                this.getDocument().remove(0, (count - maxChar));
            } catch (BadLocationException ex) {
                // Exceptions.printStackTrace(ex);
            }
        }
    }

    public void append(String s) {
        int aPos = 0;
        int aIndex = 0;
        int mIndex = 0;
        String tmpString = "";
        boolean stillSearching = true;
        String addString = remaining + s;
        remaining = "";

        if (addString.length() > 0) {
            aIndex = addString.indexOf("\u001B");
            if (aIndex == -1) {
                append(colorCurrent, addString);
                return;
            }

            if (aIndex > 0) {
                tmpString = addString.substring(0, aIndex);
                append(colorCurrent, tmpString);
                aPos = aIndex;
            }

            stillSearching = true;
            while (stillSearching) {
                mIndex = addString.indexOf("m", aPos);
                if (mIndex < 0) {
                    remaining = addString.substring(aPos, addString.length());
                    stillSearching = false;
                    continue;
                } else {
                    tmpString = addString.substring(aPos, mIndex + 1);
                    colorCurrent = getANSIColor(tmpString);
                }
                aPos = mIndex + 1;

                aIndex = addString.indexOf("\u001B", aPos);

                if (aIndex == -1) {
                    tmpString = addString.substring(aPos, addString.length());
                    append(colorCurrent, tmpString);
                    stillSearching = false;
                    continue;
                }

                tmpString = addString.substring(aPos, aIndex);
                aPos = aIndex;
                append(colorCurrent, tmpString);
            }
        }
    }

    public Color getANSIColor(String ANSIColor) {
        switch (ANSIColor) {
            case "\u001B[30m":
                return D_Black;
            case "\u001B[31m":
                return D_Red;
            case "\u001B[32m":
                return D_Green;
            case "\u001B[33m":
                return D_Yellow;
            case "\u001B[34m":
                return D_Blue;
            case "\u001B[35m":
                return D_Magenta;
            case "\u001B[36m":
                return D_Cyan;
            case "\u001B[37m":
                return D_White;
            case "\u001B[0;30m":
                return D_Black;
            case "\u001B[0;31m":
                return D_Red;
            case "\u001B[0;32m":
                return D_Green;
            case "\u001B[0;33m":
                return D_Yellow;
            case "\u001B[0;34m":
                return D_Blue;
            case "\u001B[0;35m":
                return D_Magenta;
            case "\u001B[0;36m":
                return D_Cyan;
            case "\u001B[0;37m":
                return D_White;
            case "\u001B[1;30m":
                return B_Black;
            case "\u001B[1;31m":
                return B_Red;
            case "\u001B[1;32m":
                return B_Green;
            case "\u001B[1;33m":
                return B_Yellow;
            case "\u001B[1;34m":
                return B_Blue;
            case "\u001B[1;35m":
                return B_Magenta;
            case "\u001B[1;36m":
                return B_Cyan;
            case "\u001B[1;37m":
                return B_White;
            case "\u001B[0m":
                return cReset;
            default:
                return B_White;
        }
    }

}
