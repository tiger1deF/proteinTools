Protein Tools
=======================================
Developed by Christian de Frondeville and sponsored by the Barabasi Network Science Laboratories

ProteinTools is a lightweight, flexible, and robust package that simplifies interactions with proteins. Allows for easily obtaining protein identifiers, downloading protein structural files, identifying and processing residues/residue atoms, identifying protein/ligand interactions, and much more.

ProteinTools can be downloaded via pip. ::

        pip install proteinTools

Start by creating a protein class with the desired ChEMBL, PDB, Uniprot, or HGNC/Genecard identifier (including species if not human), and use the .download method (with an optional destination directory argument) to download the PDB structural file or Alphafold representation, which generates the residues, chains, atoms, and ligands if applicable, all with their own attributes and easily accessible from the protein class.::

        from proteinTools import PT as p
    
        protein = p.Protein('1H4K')
        protein.download('/structural/protein/path')

        print(protein.Uniprot)
        print(protein.residues(44).chain)
        print(protein.residues('A_44').atoms[0].element)
        print(protein[1:3])
        print(protein.atoms(1000).center)
        print(protein.atoms(1000).line)
       
Output: ::

        P07268
        A 
        C
        [<__main__.residue object at 0x2b937f95e220>, <__main__...
        [42.103, 23.252, 48.275]
        line       ATOM   1001  OD1 ASP A 129      42.103  23.252...

.. note::

   This project is under active development.

Contents
--------

.. toctree::
      :maxdepth: 1

      installation
      documentation
