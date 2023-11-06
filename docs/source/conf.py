# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os.path
import sys

# sys.path.append(os.path.relpath("../../src"))
sys.path.insert(0, os.path.abspath("../.."))

project = 'ShockPy'
copyright = '2023, Dimitris Tsapetis'
author = 'Dimitris Tsapetis'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    "sphinxcontrib.bibtex",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery"
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

autoclass_content = "init"
add_module_names = False
autodoc_member_order = "bysource"

bibtex_bibfiles = ["bibliography.bib"]
bibtex_default_style = "unsrt"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    "logo_only": True,
    "style_nav_header_background": "#FFFFFF",
}

html_logo = "_static/logo.png"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

sphinx_gallery_conf = {
    "examples_dirs": [
        "../code",
    ],  # path to your example scripts,
    "gallery_dirs": [
        "auto_examples",
    ],  # path to where to save gallery generated output
    "binder": {
        # Required keys
        "org": "dimtsap",
        "repo": "ImpedancePy",
        "branch": "main",  # Can be any branch, tag, or commit hash. Use a branch that hosts your docs.
        "binderhub_url": "https://mybinder.org",
        # Any URL of a binderhub deployment. Must be full URL (e.g. https://mybinder.org).
        "dependencies": "../../requirements.txt",
        "notebooks_dir": "notebooks",
        "use_jupyter_lab": True
        # Jupyter notebooks for Binder will be copied to this directory (relative to built documentation root).
    },
    "ignore_pattern": "/local_",
}