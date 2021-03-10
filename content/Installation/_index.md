+++
title = "Installation"
date = 2021-02-23T10:07:43+05:30
weight = 5
chapter = false
+++

# DJMol

## Download and Installation
---
 
### Download
+ #### [Binaries](https://djmol.s3.amazonaws.com/DJMolPlatformv2.1_Win64.zip)  
+ #### [Source](https://djmol.s3.amazonaws.com/DJMolSource_2020_BETAformat.zip)

## For Compiling/Installation of DJMol and Add-Ons

### General Information

A familiarity with programming using Netbeans IDE is assumed. All of the DJMol Program was compiled by using Netbeans IDE, v 8.2 (64bit) in  Windows 64 OS.

The needed Packages are:
:---
+ JAVA 1.8 (jdk-8u201-windows-x64.exe)
+ Netbeans 8.2 (netbeans-8.2-javaee-windows.exe)
+ Python plug-in (2017-08-29-nbpython-nbms.zip)
+ Python 3.7.2 (python-3.7.2-amd64.exe)
+ Pip script (get-pip.py)
+ Numpy (numpy-1.16.0+mkl-cp37-cp37m-win_amd64.whl)


We also used (64 bit) the following Python packages with version numbers:
:---
+ SciPy (scipy 1.2.1)  
+ ASE (ASE 3.17.0)
+ Matplotlib (Matplotlib 3.0.3)

> And these modules can be downloaded by using pip (use pip3).


A. Compiling Main Program
:---

1. Open the Project Directory, ***DJmol Platform v2.1*** and Compile and Build DJMol program. See the YouTube channel for its demonstration.
2. Then make a Zip Distribution (by Selecting project ***DJmol Platform v2.1*** followed by selecting  ***Package as***, ***Zip Distribution***). A new ***dist*** folder will be created and it holds the Zip file.

B. Installing Program
:---

3. Go to ***dist*** folder and unzip the file and Copy all the files/folders from the ***Auxiliary*** folder into *** dist\djmolplatform1***. This is shown in the demonstration.
4. Go to ***dist\djmolplatform1\etc*** and replace the following line in the ***djmolplatform1.CONF*** file:
        `default_options="--branding djmolplatform1 -J-Xms24m -J-Xmx64m" `   to
        `default_options="--branding djmolplatform1 -J-Xms240m -J-Xmx640m"`
5. To adjust the Resolution of the program  Right click on  ***djmolplatform164.exe*** then move to   ***Compatibility*** tab and click on *** Change high DPI settings*** and Select the option "Override high DPI scaling behavior".
6. Follow the instruction in the PathVariablesSetting.pdf file.

DJMol program can be now simply executed by double clicking ***djmolplatform164.exe*** which is located in ***bin*** folder.  

Note that for the First Time of the execution, you may want to select,  ***Disable Modules and Continue*** option. To fix this, you can try to clean your user-directory as it is mentioned in http://wiki.netbeans.org/FaqWhatIsUserdir

And once the application starts do the following (not mandate but it is useful):

```
Close <Start Page>
Close <Services>
Add <Windows -> Favorites>
Add <Windows -> Output>
```

C. Compiling Add-On Programs
:---

See ***AddOns_compilation.txt*** in Add-On directory.
