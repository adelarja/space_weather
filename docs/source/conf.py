# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import datetime
import os
import solarwindpy
import sys

# this path is pointing to project

sys.path.insert(0, str(solarwindpy))
#sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------

project = 'solarwindpy'
copyright = '2022, Adriana Gulisano, Adel Arja, Violeta Bazzano, Ricardo Pafundi'
author = 'Adriana Gulisano, Adel Arja, Violeta Bazzano, Ricardo Pafundi'

# The full version, including alpha/beta/rc tags
release = solarwindpy.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    'sphinx.ext.autodoc', #permite crear doc automaticamente de docstrings
    'sphinx.ext.intersphinx',# crea enlaces a documentos
    'sphinx.ext.coverage',
    'sphinx.ext.todo',# ToDo cuentas por hacer
    'sphinx.ext.mathjax',# Permite agregar formulas mat escritas en LaTeX
    'nbsphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode', #codigo se incluye en la documentacion
    'sphinx_rtd_theme'     #permite usar plantillas html
]

#autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.

#autodoc_mock_imports = ['bs4','requests'] #importa el codigo de ejemplo

# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
#language = 'y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
#html_sidebars = { '**': ['globaltoc.html','relations.html', #Contenidos Barra Lateral
#            'sourcelink.html','searchbox.html'], }
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

master_doc = 'index'

