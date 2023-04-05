Protein Tools
=======================================
Developed by Christian de Frondeville and sponsored by the Barabasi Network Science Laboratories

.. image:: https://github.com/ChatterjeeAyan/AI-Bind/blob/main/Images/NetSci_Logo.png
   :width: 200px
   :height: 100px
   :scale: 50 %
   :alt: alternate text
   :align: right

ProteinTools is a lightweight, flexible, and robust package that simplifies interactions with proteins. Allows for easily obtaining protein identifiers, downloading protein structural files, identifying and processing residues/residue atoms, identifying protein/ligand interactions, and much more.

ProteinTools can be installed via the following command: ::

        pip install proteinTools
        
Start by creating a protein class with the desired ChEMBL, PDB, Uniprot, or HGNC/Genecard identifier (including species if not human), and use the .download method (with an optional destination directory argument) to download the PDB structural file or Alphafold representation, which generates the residues, chains, atoms, and ligands if applicable, all with their own attributes and easily accessible from the protein class.::

        import proteinTools
        
        protein = Protein('1H4K')
        print(protein.Uniprot)
        print(protein.residues('A44')
        print(protein[1:3])
        print(protein.atoms(1000))
        print(protein.FASTA)
       
Output: ::

        P07268
        Chain                                                         A
        Amino Acid                                                  ASN
        Index                                                        44
        Atoms         [<__main__.atom object at 0x2b937fb49eb0>, <__...
        [<__main__.residue object at 0x2b937f95e220>, <__main__.residue object at 0x2b937f95ee80>]
        Element                                                    O
       x                                                     42.103
       y                                                     23.252
       z                                                     48.275
       Residue                                                 A129
       Line       ATOM   1001  OD1 ASP A 129      42.103  23.252...      

HKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAVARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALG



Full documentation is available at
https://proteintools.readthedocs.io/en/latest/
