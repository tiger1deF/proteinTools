Usage
=====

.. _installation:

Installation
------------

To use proteinTools, first install it using pip:

.. code-block:: bash

   (.venv) $ pip install proteinTools

Creating Proteins
----------------

To create a protein item, you can use the ``Protein(<Identifier>, <species>)`` method, where 
species is human by default. Identifiers can be ChEMBL IDs, UniprotKB/Uniprot IDs, PDB IDs, or HGNC/Genecard IDs. 

The protein file can than be downloaded to a directory (default the user's current directory) with the ``.download`` method. Upon downloading the protein file (.pdb/.mmCIF from RCSB in the case of a PDB ID, or an Alphafold-generated structure otherwise), the protein will be populated with chains, residues, atoms, and other identifiers from the file.

.. code-block:: python

   import proteinTools as p
   myprot = p.Protein('2D3Z', 'e coli')
   myprot.download('/path/to/directory')
   
Protein Properties
------------------
The protein's ChEMBL ID, Uniprot ID, PDB ID, and HGNC/Gene ID can all be accessed with properties of the same name

.. code-block:: python
   
   myprot.Uniprot
   myprot.HGNC
   myprot.ChEMBL

Residues can be queried by indexing, or by using the residues method.

.. code-block:: python

   myprot[1:3]
   myprot[100]
   myprot.residues('A_55')
   
Protein interactions can be accessed with the property of the same name, returning BindingDB ligands and their activities, ChEMBL ligands and their activities (all in uM) and the ID of STITCH ligands and proteins, all returned in a dictionary.

.. code-block:: python
   
   myprot.interactions
   
The total amount of residues in the protein is obtainable simply by using the len() magic method.

.. code-block:: python

   protein_length = len(myprot)
   
Residue Properties
-------------------

Residue amino acids (AA), chain, atoms, index, and name can be accessed by properties of the same title.

.. code-block:: python

   protein[1].name
   protein.residues('A433')['Name']
   protein[5].AA
   protein[8].chain
   protein[2].atoms

The center of mass of each residue can be calculated with the ``.center`` property.

Atom Properties
----------------

