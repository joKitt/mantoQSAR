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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.io.IOException;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;


final class ConsolePanel extends JPanel implements Runnable{

    private ColorTextPane textArea;
    private final Thread reader;
    private final Thread reader2;
    private final boolean quit; 
    
    private final PipedInputStream pin=new PipedInputStream();
    private final PipedInputStream pin2=new PipedInputStream();

    public ConsolePanel(){

        
        textArea=new ColorTextPane();
        textArea.setMargin(new Insets(10,10,10,10));
        
        this.setBackground(Color.WHITE);
        this.setLayout(new BorderLayout());
        this.add(new JScrollPane(textArea), BorderLayout.CENTER);
 
        this.appendMsg("mantoQSAR version 1.0 alpha \n");
        
        try{
            PipedOutputStream pout=new PipedOutputStream(this.pin);
            System.setOut(new PrintStream(pout,true));
	}catch (java.io.IOException io){
            this.appendMsg("Couldn't redirect STDOUT to this console\n"+io.getMessage());
	}catch (SecurityException se){
            this.appendMsg("Couldn't redirect STDOUT to this console\n"+se.getMessage());
        }

        try{
            PipedOutputStream pout2=new PipedOutputStream(this.pin2);
            System.setErr(new PrintStream(pout2,true));
	}catch (java.io.IOException io){
                this.appendMsg("Couldn't redirect STDERR to this console\n"+io.getMessage());
	}catch (SecurityException se){
                this.appendMsg("Couldn't redirect STDERR to this console\n"+se.getMessage());
		}

		quit=false; 

		reader=new Thread(this);
		reader.setDaemon(true);
		reader.start();

		reader2=new Thread(this);
		reader2.setDaemon(true);
		reader2.start();
    
}

    public synchronized void actionPerformed(ActionEvent evt)
	{
		textArea.setText("");
	}
    

    @Override
    public synchronized void run() {
       
        try
		{
			while (Thread.currentThread()==reader)
			{
				try { this.wait(100);}catch(InterruptedException ie) {}
				if (pin.available()!=0)
				{
					final String input=this.readLine(pin);
                                        this.appendMsg(input);
					
				}
				if (quit) return;
			}

			while (Thread.currentThread()==reader2)
			{
				try { this.wait(100);}catch(InterruptedException ie) {}
				if (pin2.available()!=0)
				{
					String input=this.readLine(pin2);
                                        this.appendMsg(input);
					// textArea.append(input);
				}
				if (quit){
                                    return;
                                } 
                                    
			}
		} catch (Exception e)
		{
                    this.appendMsg("\nConsole reports an Internal error.");
                    this.appendMsg("The error is: "+e);
		}

    }

    
    public synchronized String readLine(PipedInputStream in) throws IOException
	{
		String input="";
		do
		{
			int available=in.available();
			if (available==0) break;
			byte b[]=new byte[available];
		in.read(b);
			input=input+new String(b,0,b.length);
		}while( !input.endsWith("\n") &&  !input.endsWith("\r\n") && !quit);
		return input;
	}

    public void appendMsg(final String str) {

        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                textArea.append(str);
            }
        });
    }
   }
