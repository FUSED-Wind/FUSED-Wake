from setuptools import setup, find_packages

#####################################################################
kwargs = {'author': 'Juan Pablo Murcia, Pierre-Elouan Rethore, DTU Wind Energy',
        'author_email': 'jumu@dtu.dk',
        'classifiers': ['Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering'],
        'description': 'Framework for Unified Systems Engineering and Design of Wind Plants',
        'include_package_data': True,
        'keywords': ['openmdao'],
        'license': 'Apache v. 2.0',
        #'maintainer': '',
        #'maintainer_email': '',
        'name': 'fusedwake',
        'packages': find_packages('*'),
        #'url': '',
        'version': '0.1',
         'install_requires': ['openmdao.main', 'numpy', 'scipy', 'pandas', 'matplotlib', 'seaborn','fusedwind'],
         'dependency_links': ['git+http://github.com/FUSED-Wind/fusedwind.git@0.1.0#egg=fusedwind'],
      'zip_safe': False}

setup(**kwargs)

##################################################################


