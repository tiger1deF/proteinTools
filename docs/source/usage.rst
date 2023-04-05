Usage
=====

.. _installation:

Installation
------------

To use Lumache, first install it using pip:

.. code-block:: console

   (.venv) $ pip install proteinTools

Creating Proteins
----------------

To create a protein item, you can use the ```Protein(<Identifier>, <species>) method, where 
species is human by default:

.. autofunction:: proteinTools.Protein('2D3Z', 'e coli')


.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

