+++
title = "Chaper 6"
weight = 5
chapter = true
+++

# Tools of  DJMol

Tools are generally calculator independent (add-on) modules which can
assist some steps in the modeling.

Some of the important Tools are described below.


### [1] Analysing Point Group

DJMol use a point group analyzer and from this one can easily find the
point group of a molecule and its subgroups. It can also tell the next
possible point group within a specified error-bar. This module can be
invoked by:

Tools ➤ PointGroup Analyzer

This tool can also be used to deform the geometry of molecules so that
one can modify a molecule’s PG into another one (in other words, we can
constrain a PG into its sub group).

![ ](/images/image33.png " ")

### [2] Z-Matrix Editor

           Z-Matrix or Z-mat tool can be called by:

Extra ➤ Z-Matrix Tool

This Z-Mat editor is made for the XYZ files; It can be either used as an
independent program (use File button to load a required XYZ file) or as
an auxiliary program to manipulate the loaded XYZ file. Use either \<New
File\> button to vary Z-matrix a of a new Cartesian file or \<Current
File\> button to fetch the currently loaded geometry from the main panel
of DJmol.

![ ](/images/image34.png " ")

Figure 6.1: A butane molecule (exaggerated) and its Z-matrix data.

To change Z-Mat Data:

After displaying the geometry, select \<Select a Z-mat variable\> radio
button to show the corresponding Z-matrix data (MOPAC style, generated
from the OpenBabel application) of the displayed molecule in this
application.

 

Also select any single cell from either \<X\>, or \<Y\> or \<Z\> column,
which represents bond-length, bond angles, or dihedral angles,
respectively. After that select \<Insert Variable\> button to make the
cell element as a Z-mat variable.

And then use the \<Save\> button and close the window. After that,
immediately select the \<Start Zmatrix Variation\> menu from the top of
the window, and by using [+] and [-] buttons, start the variation.

One can either save the data (\<Update Zmat\> button) or discard the
variation and to move the molecule into its original state by using
\<Undo Update\> button.

It is also important to note that to stop and undo the current
variation,first to select \<Stop Variations\> and press \<Undo Update\>
button\>.

If you want to move only X or Y or Z positions of a single atom use the
menu ( \<Move selected Atom\> ), and select an atom (use a click on the
desired atom) , and press any[x+][x-]...buttons near the bottom left
corner.

Important: To save the updated Geometry press \<Update Zmat\> AND then
\<Save Structure\> button.

Or to discard the variation and to move the molecule into its original
state by using \<Undo Update\> button.

Also don't forget to use \<Save Structure\> button to save the modified
Cartesian file. It will be saved as CurrentFinalXYZ.xyz file in the
folder \< ZMATfolder\> in the program directory.

If necessary use:

[ 1 ] \<Select Atoms Labels\>  to see individual atom labels, which may
be useful when one manipulates Z-Mat.

[ 2 ] \<Show initial bonds\> for keeping the original bonding pattern,
irrespective of the current geometric state. This can be useful when
some bonds are too far away from the initial bonded state.

[ 3 ] \<Show axis\> to show cartesian axis

[ 4 ] \<Add New Bond\>[n1] [n2] - to show a new bond between the
atomlabel, n1 to atomlabel2, n2 atoms.

[ 5 ] \<Builder\> option can be used if necessary.

Note that a user can also enter numerical values (in float) to the Text
Fields. Also Move the slide bar to maximum when you vary the
Angles/Dihedral angles and to a minimum if you want to manipulate the
distance.

### [3] Remote Submission of Scripts

For time consuming job, such as MD calculations, it is a better idea to
use a dedicated machine like Linux cluster/Google Cloud etc. instead of
a PC. In this case, one need to generate SSH keys (ie. one private and
one public key), preferably, by using PUTTYGEN application; see
[Puttygen](https://www.puttygen.com/) to download this mandate application.

Note that we need \<RSA\> based private key (which is stored in your PC)
and a public key (stored in the Remote machine). Note that you must
enter a \<Passphrase\> for the \<Private key\> and its format should be
in the OpenSSH (ie. not in the PPK format, which is the PuttyGens'
default). PuttyGEN can convert this private key in PPK format to OpenSSH
format. Finally the Private Key should be kept in the \<ssh\> folder
with a file name, \<id\_rsa\>. Make sure that the file \<known\_hosts\>
is existed in \<ssh\> format. After this one can invoke the \<Remote
Submission Tool\> as:

RunProject ➤ Remote Submission

This tool can be used instead of PUTTY, and it is equipped with tools
for file upload/download purposes. The downloaded files are kept in
\<SCP\> folder, by default.

        ![ ](/images/image35.png " ")

Figure 6.2: A sample SSH session (for Google Cloud Platform remote
connection and file transactions).

### [4] Version Control using GIT {.c75 .c13}

Version control (VC) systems are a category of software tools that help
a software team manage changes to source code over time. Version control
software keeps track of every modification to the code in a special kind
of database. If a mistake is made, developers can turn back the clock
and compare earlier versions of the code to help fix the mistake while
minimizing disruption to all team members. Git is an example of VCS.

Why the VC is used in the DJMol? Since the text based input scripts are
playing an important role in the modeling (inputs of Siesta/DFTB/ASE,
all are text files). Using VC, these data can be re-edited remotely and
at the same time it keeps different version of these files in a
systematic way (for example, later, it can be used by
collaborators/public).

Some Basic Definitions to be familiarized are:

​1.  FETCH: The git fetch command downloads commits, files, and refs
from a remote repository into your local repo. Fetching is what you do
when you want to see what everybody else has been working on.

​2.  PULL: The git pull command is used to fetch and download content
from a remote repository and immediately update the local repository to
match that content.

​3.  PUSH: The git push command is used to upload local repository
content to a        remote repository. Pushing is how you transfer
commits from your local repository to a remote repo. It is the
counterparts to git fetch.

How to use Git to push and pull in DJMol?

​1.  Initializing Version Control.

​a.  Set up a Git Remote Repository in GitHub and copy the Repository
URL.

​b.  Open the project in which you need to use version Control.

​c.  Initializing the Local Repository.

Team ➤ Git ➤ Initialize ➤ (OK)

​d.  Linking the Local Repository and the online Repository

Team ➤ Remote ➤ Pull ➤ Fill ➤ Check Master ➤ Finish

Specify  ➤ PateTheRepository ➤ UserName/PassWord (if needed) ➤  Next

         2.  Push

​a.  After making necessary changes in the Project.

​b.  From the Project sidebar, add the needed files to the repository by
right-clicking on the file and Clicking on Add.

​c.  Commit the files to the Repository

Team ➤ Commit ➤ AddACommitMessage ➤  Commit

​d.  Push to remote repository.

Team ➤ Remote ➤ PushToUpstream ➤  Yes

​3.  Pull

​a.  Fetch

Team ➤ Remote ➤ Fetch ➤ Next ➤ Finish

​b.  Pull

Team ➤ Remote ➤ Pull ➤ Next ➤ CheckMaster ➤ Finish

Installation of Git

Adding Github Plugins

After installing the application, search \<Git\> in the plug-ins search
field. Select the \<Git\> from the \<Installed Packages\> and click on
\<Activate\> as shown below figure.

![ ](/images/image30.png " ")

Figure 6.3: Installation of Github plugin in DJMol.


### [5] Database File Retrieving {.c75 .c13}

To fetch different structural file form various open file repositories,
a Database tool is implemented. It is essentially uses Jmol’s DB module.

It can be launched by:

Extra ➤ Database Tool

Open URL: It can be used to open XYZ or PDB file or any other Jmol
recognizable file from the web server, say from Github.  

Open MOL: It is mainly for, SMILES, InChI, or CAS from either a PubChem
database or from NCI/NIH database.

Open PDB: It use RCSB web (please prefer RCSB) to load 4-character PDB
ID (eg. 1crn) or 3 -character ligand (eg. 60C).

Open COD: It opens a specific COD ID from
http://www.crystallography.net/cod/ 

Open Materials Project: It opens a specific Materials Project ID number:

All opened file is saved in ./Database folder from  

Export To: Only MOL, XYZ and PDB formats are supported for export its
images.

![ ](/images/image21.png " ")

Figure 6.4: Showing a Database window with a retrieved COD file from COD
databse.

### [6] Miscellaneous Tools {.c75 .c13}

Some of the miscellaneous tools are shown below:

Start Batch Process: For DFTB+ calculation, a number of HSD files in a
particular folder can be called one after another (known as batch 
processing). This will be useful for building 1D PES or to analyze
energies of particular set of molecules. The resultant TGZ files contain
all the out files including the submitted HSD file, and it will be
transformed to the folder at the end of the each calculation.

![ ](/images/image22.png " ")

Calculator: System’s default calculator can be invoked.

GIF Utility: GIF animation files can be viewed or generated using this
utility (eg. phonon modes). Time between two frames can be adjusted and
the resultant file is stored in ./Scratch\_images.

Process Status: It is used to monitor current resources of the computer
system (used/availble RAM, HDD space etc.). This can be used before
starting a calculation. For example, if the PC memory is low, allocate
more RAM by stopping other less important process by using Window’s Task
Manager utility.

Unit Conversion: A basic unit converter for length, time, energy etc. it
also includes atomic unit.

Matrix Viewer: 2D matrix data viewer of the program. It is useful to
analyze the overall shape and symmetry of the matrix data, such as
Hessian or Hamiltonian. See example directory for its sample file.

Convert Structure: Open Babel is used to convert from one structure data
into another data. See Figure 2.10.

Standard Orientation: Using symmetry a disoriented molecule can be
oriented with respect to a symmetry axis. This tool is useful for making
a more systematic Z-matrix or Cartesian file . It is strongly
recommended that the standard oriented geometry data should be used in
the Z-matrix tool to re-adjust its coordinates.  

![ ](/images/image23.png " ")

Open DJMol Directory: It will open the parent directory.

Send Message: User can mail the forum using this utility. Please
academic E-mail ID (avoid .com E-mail IDs – it won’t support). If needed
images should be linked as an URL link (say, by drive.google.com/…
link).
