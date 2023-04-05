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

   myprot = proteinTools.Protein('2D3Z', 'e coli')
   myprot.download('/path/to/directory')
   
.. autofunction:: to_csv

