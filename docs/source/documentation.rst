================
Full Documentation
================

.. note:: Python 3 is required to run proteinTools. Any version of python 2 is currently not supported.

Creating Proteins
----------------

To create a protein structure, you can use the ``Protein(<Identifier>, <protein-ID> <species>)`` initializer, where 
species is optional (human by default) and protein-ID is also optional (PDB by default). Identifiers can be ChEMBL IDs, UniprotKB/Uniprot IDs, PDB IDs, or HGNC/Genecard IDs.

The protein file can than be downloaded to a directory (default the user's current directory) with the ``.download`` method. Upon downloading the protein file (.pdb/.mmCIF from RCSB in the case of a PDB ID, or an Alphafold-generated structure otherwise), the protein will be populated with chains, residues, atoms, and other identifiers from the file.

.. code-block:: python

   from proteinTools import PT as p
   myprot = p.Protein('2D3Z', 'e coli')
   myprot.download('/path/to/directory')
   
   newprot = p.Protein('P02786', 'Uniprot', 'human')
   myprot.download('/path/to/directory')
   
   newprot = p.Protein('L2HGDH', 'HGNC')
   myprot.download('/path/to/directory')
   
Protein Properties
------------------
The protein's ChEMBL ID, Uniprot ID, PDB ID, and Gene ID can all be accessed with properties of the same name. For PDB IDs, proteinTools returns a **dataframe** containing the PDB ID(s), minimum resolution of the PDB file in Angstroms, and the total number of unique ligands present in the file. The protein file does not need to be downloaded to access these values.

.. code-block:: python
   
   myprot.Uniprot
   myprot.Gene
   myprot.ChEMBL
   myprot.PDB

Protein residues can be queried by indexing as one would a list.

.. code-block:: python

   myprot[100:350]
   myprot[100]
   
Residues and atoms can also be accessed via the .atoms and .residues method, which return the corrosponding atom and residue object.

.. code-block:: python
   
   myprot.residues('A_144')
   myprot.atom(6572)
   
The protein's FASTA sequence can be easily accessed by the FASTA property, which generates a FASTA sequence directly from the file.

.. code-block:: python

   myprot.FASTA
   
Protein-ligand interactions can be accessed with the ligand_interactions property, returning a dictionary containing BindingDB ligands and their activities, ChEMBL ligands and their activities (all in uM) and the ID of STITCH ligands and proteins. The protein does not need to be downloaded to access these values.

.. code-block:: python
   
   myprot.ligand_interactions

Protein-protein interactions can be accessed with the protein_interactions property, returning a dictionary containing stringDB protein keys with a list containing the maximum interaction score for the PPI and mean interaction score for the PPI. The protein does not need to be downloaded to access these values.

.. code-block:: python

    myprot.protein_interactions

The total amount of residues in the protein is obtainable simply by using the len() magic method.

.. code-block:: python

   protein_length = len(myprot)
   
A list of every atom in the protein and their properties can be created with the .to_csv(<destination>) method, where the default destination is the user's current directory.

.. code-block:: python
  
   myprot.to_csv('/path/to/directory')

Every line of the file is present in the .data method.

.. code-block:: python

   protein_lines = myprot.data


the .ligands property returns a two-row dataframe, with the first row containing the primary ligand, all unique cofactors, and all ions present in the protein structure, where the second row contains the corrosponding ligand objects populated directly from the file for each category. 

.. code-block:: python

   unique_ligands = myprot.ligands

Ligands can be queried from the .ligand_list attribute, which returns a list of every unique ligand position of every ligand present in the structural file (outside of water moelecules). Ligand properties are described below.

.. code-block:: python

   for ligand in myprot.ligand_list:
         ligands.append(ligand)
         
.. Note:: 

   CIF (mmCIF) files are currently supported by proteinTools, but certain functionalities (stripping ligand sites and secondary structure) are not available.
   
Residue Properties
-------------------

Residue amino acids (AA), chain, atoms, index, and name can be accessed by properties of the same title.

.. code-block:: python

   myprot[1].name
   myprot.residues('A433')['name']
   myprot[5].AA
   myprot[8].chain
   myprot[2].atoms
   residues = myprot[1:100]

The center of mass of each residue can be calculated with the ``.center`` property, which returns a list of the x, y, and z coordinate of the residue center.

.. code-block:: python

   residue_center = myprot[1].center
   
If the protein has a PDB ID format, the secondary structure of each residue can also be obtained with the structure property (between HELIX, SHEET, and UNSTRUCTURED). 

.. code-block:: python

   residue_structure = myprot[160].structure

Atom Properties
----------------

The x, y, and z coordinate of atoms, as well as their mass, center element, line (line data from protein file), and atomic radius can be accessed by properties of the same title.

.. code-block:: python

   residue, elements, positions = myprot.residue('B123'), [], []
   for atom in residue.atoms:
        elements.append(atom.element)
        atom_volume = 4 / 3 * math.pi * pow(atom.radius, 3)
        positions.append(atom.center)
        
The parent residue for atoms constructed via a protein file can be accessed with the ``.parent_residues`` property.

.. code-block:: python
   
   for atoms in residue.atoms:
      protein_chain = atoms.parent_residue.chain  
    
Ligand Properties
--------------
If the protein is a PDB file containing ligands (that are not water molecules), they will automatically be added to the .ligands protein attribute. The ligand ID as present in the PDB file can be accessed with the ID attribute, and atoms of the atom class can be accessed with the atoms attribute.

The center of mass and the radius of gyration of each ligand can be calculated via their respective properties.

.. code-block:: python

   ligand.center
   ligand.radius

The ligand file can be downloaded by the ``.download('/path/to/file')`` method, which defaults to the user's current directory and saves the ligand in .sdf format.

.. code-block:: python

   ligand = protein.ligand_list[3]
   ligand.download()
   for ligand in protein.ligand_list:
       print(ligand.ID)
       print(ligand.center)

Ligand objects can also be created by specifying the path of the existing file when instantializing the ligand. Acceptable file formats are .pdb, .sdf, or .pdbqt. 

.. code-block:: python

      ligand = protein.ligand(file_path = '/path/to/file')
      center = ligand.center

Ligand files can also be instantiated separate of the protein. Simply generate the ligand with the ligand ID, and use the .download method with the path to download (defaults to the current user directory). InChiKeys, PDB IDs, and SMILEs sequences are accepted.

.. code-block:: python

    lig = p.ligand('C1=CC=C2C(=C1)C=CC=C2CCC(CO)N3C=C(N=C3)C(=O)N')
    lig.download('/path/to/directory')
    
    
.. warning:: 

   The above functionality depends on the chemical structure package Openbabel, which can be installed via ``pip install openbabel``

Autodock Vina Integration
--------------
The protein file can be downloaded in .pdbqt format (a specialized Vina protein file format for docking) using the .pdbqt('/path/to/directory') method, where the path defaults to the current user directory.

.. code-block:: python

      protein.pdbqt()

A recreation of the protein file stripped of ligand molecules can be created using the .strip_ligands method, where the destination is the location of the file, defaulting to the current user directory.

.. code-block:: python

   myprot.strip_ligands()

A recreation of the protein file stripped of only the primary ligand molecules (not cofactors and ions) can be created using the .strip_primary method, where the destination is the location of the file, defaulting to the current user directory.

.. code-block:: python

   myprot.strip_primary()
