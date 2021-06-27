.. title:: Installation of DJMOL

===========================================
Installation of DJMOL
===========================================

Download 
===========================

        - `Binaries V2.1 (June 2021) <https://drive.google.com/file/d/1whpwssDTE43mL4HGfY4CsSo4IRKqBC_S/view?usp=sharing>`_
        - `Source V2.1 (June 2021) <https://drive.google.com/file/d/1jRCYA33f3bsJussPG9fqxCGmxMlocsmR/view?usp=sharing>`_
  


Installation of Binary
=======================

General Information:
-----------------------

This descibes how to install Binary of DJMol.

The following programs are required (Mandate) to Run the software
    
- Java 1.8 or higher (if it is not installed please install
  Java via installing JDK from 
  https://www.oracle.com/technetwork/java/javase/downloads/index.html)
  or from http://www.djmol.info/download.html 
- Python 3.X (3.7 preferred)
  Please use 64 bit version of Python. 
- NumPy, ASE and (optionally) Python3 plugg-in for Netbeans.


Running DJMol v 2.1 in Windows 10 OS 64 bit
---------------------------------------------

- Unzip the file, <dmolplatform1.zip> and Copy all the files/folders from the
  <./Auxilary> folder into <./djmolplatform1/.>. 
  See https://www.youtube.com/watch?v=zm6Mbh0T1m8 for the demonstration.

- To adjust the Resolution of the program  Right click on  <djmolplatform164.exe> then move to 
  <Compatibility> tab and click on <Change high DPI settings> and Select the option "Override high DPI scaling behavior".

- Please follow the instruction in the <PathVariablesSetting_v2.pdf> file.

DJMol program can be now simply executed by double clicking <djmolplatform164.exe>
which is located in <bin> folder.  

.. note:: 

    Note that for the First Time of the execution, you may want to select, 
    <Disable Modules and Continue> option.
    To fix this, you can try to clean your user-directory as it is mentioned in 
    http://wiki.netbeans.org/FaqWhatIsUserdir 
        

Compilation from the source code
=====================================


General-Information
-----------------------

All of the DJMol Program was compiled by using Netbeans IDE, v 8.2 (64bit) for Windows 64 OS.

This program and all other mandate packages for compiling the software are availble at free of cost and
to download these packages please see: https://djmolprogram.github.io

-	JAVA1.8		jdk-8u201-windows-x64.exe
-	Netbeans 8.2	netbeans-8.2-javaee-windows.exe
-	Python pluggin	2017-08-29-nbpython-nbms.zip
-	Python 3.7.2	python-3.7.2-amd64.exe
-	Pip script 	get-pip.py
-	Numpy		numpy-1.16.0+mkl-cp37-cp37m-win_amd64.whl


We also used (64 bit) python packages:
-	scipy 1.2.1  
-	numpy 1.16.0 
-	ASE 3.17.0
-	Matplotlib 3.0.3 

And these programs can be downloaded by using pip.


Compiling Main-Program <DJMolplatform>
----------------------------------------

- Open the Project Directory, <DJMolCompPhysComm_2020_BETAformat> and Compile and Build DJMol program.
    See https://www.youtube.com/watch?v=9gQzd9qnjN0 for its demonstration.

- Then make a Zip Distribution (by Selecting project <DJmol Platform v2.1> followed by selecting 
    <Package as>, <Zip Distribution>). A new <dist> folder will be created and it holds the Zip file.

- Go to <dist> folder and unzip the file and Copy all the files/folders from the
  <Auxilary> folder into <dist\djmolplatform1>. 
  See https://www.youtube.com/watch?v=zm6Mbh0T1m8 for the demonstration.

- Go to <dist\djmolplatform1\etc> and replace the the following line in the <djmolplatform1.CONF> file:
  
  default_options="--branding djmolplatform1 -J-Xms24m -J-Xmx64m"
    
  to
	
  default_options="--branding djmolplatform1 -J-Xms240m -J-Xmx640m" (meaning RAM is defined between 240-640MB for DJMol)

- To adjust the Resolution of the program  Right click on  <djmolplatform164.exe> then move to 
    <Compatibility> tab and click on <Change high DPI settings> and Select the option "Override high DPI scaling behavior".

- Please follow the instruction in the PathVariablesSetting_v2.pdf file.

DJMol program can be now simply executed by double clicking <djmolplatform164.exe>
which is located in <bin> folder.  

Note that for the First Time of the execution, you may want to select, 
<Disable Modules and Continue> option.
To fix this, you can try to clean your user-directory as it is mentioned in 
http://wiki.netbeans.org/FaqWhatIsUserdir


Compiling Add-On Programs:
-----------------------------
See <AddOns_compilation.txt> in Add-On directory.

For Furthur info: Please see the <DJMol Channel> on <Compilation> and <Installation>.
    
System Requirements
=======================

* Software:
    * Windows 7/10 64 bit OS
    * Java 1.8 or later, Python 3.x

* Hardware:
    * At least 20 GB free space to install DJMol and add-ons.
    * Preferred RAM size is 4 GB or higher.
