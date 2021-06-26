.. title:: Chapter-5 :: HTML documentation

======================================
Chapter 5: Python scripting with ASE
======================================

To demonstrate ASE scripting within the platform, we must create a Python project:

**File ➤ New Project ➤ Python ➤ Python Project – Ant** 

and save this in a suitable directory. We name the main class as <ASEdemo>.

.. figure:: _static/webdocimages/036.png
   :alt: 036
   :align: center
   :scale: 80%

   Setting up a Python ASE project in the DJMol application.

.. figure:: _static/webdocimages/037.png
   :alt: 037
   :align: center
   :scale: 70%

   A script window displays an ASE script.

The purpose of this ASE calculation is to obtain the phonon dispersion for bulk aluminum using 
a 7×7×7 supercell within effective medium theory (EMT as it is implemented in ASE). For 
the source, see, Ex2_phononbulkAl_DOS.py file from the example directory. After running it, 
this script will create the band and the DOS diagram of the system. 

The script is provided in the Example directory. Once the script is ready, use <Run Project> 
command from the <Run> menu to execute the script. It will initiate ASE engine and produce 
the band and DOS diagrams. Optionally, users can develop their own python scripts in this 
project directory for other tasks (eg. to save the animation of a particular mode).

**Run ➤ Run Project**

.. figure:: _static/webdocimages/038.png
   :alt: 038
   :align: center
   :scale: 70%
   
   (1) is the Workplace and (2) is the console output from the Python run and 
   (3) is the resultant Band and DOS figures that is created in the Workplace.

In the next advanced example (see example directory of Python) we will show that RDF 
(radial distribution function) of melting copper in a cubic box. Here, EMT is also used but RDF 
function is taking from ASAP calculator. ASAP is a calculator for doing large-scale classical 
molecular dynamics within the ASE package (its libraries can be compiled in CygWin or WSL). 
Below figure shows the RDF from ASAP calculator and uses systems terminal.

.. figure:: _static/webdocimages/039.png
   :alt: 039
   :align: center
   :scale: 70%
   
   Radial distribution function of melting copper (FCC) at 2500K from an ASAP calculator.

In first add Terminal into the system by:

Window ➤ IDE Tools ➤ Terminal 

This will get CygWin ot WSL, depends upon the settings. After this, type the code in the terminal:

.. code:: python
    
    python3.6m asedemo.py

It will launch an ASAP calculation, and calculates RDF of the melting copper atoms. 
Note that its RDF is very similar to that of water and it confirms that liquid state of Cu 
atoms at high temperature is well represented by the ASAP model.