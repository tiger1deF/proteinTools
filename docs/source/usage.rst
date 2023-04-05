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
species is human by default:

.. code-block:: console

   proteinTools.Protein('2D3Z', 'e coli')


