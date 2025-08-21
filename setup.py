#region modules
from setuptools import setup, find_packages
#endregion

#region variables
#endregion

#region functions
setup(
    name='fpflow',
    version='0.1.0',
    description='First principles workflow',
    long_description='First principles workflow',
    author='Krishnaa Vadivel',
    author_email='krishnaa.vadivel@yale.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={'fp': ['data/**/*']},
    requires=[
        'numpy',
        'scipy',
        'ase',
        'dill',
        'pyyaml',
    ],
)
#endregion

#region classes
#endregion