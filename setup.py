#region modules
from setuptools import setup, find_packages
#endregion

#region variables
#endregion

#region functions
setup(
    name='fpflow',
    version='1.1.0',
    description='First principles workflow',
    long_description='First principles workflow',
    author='Krishnaa Vadivel',
    author_email='krishnaa.vadivel@yale.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    package_data={'fpflow': ['data/**/*']},
    install_requires=[
        'numpy',
        'scipy',
        'ase',
        'pyyaml',
        'jmespath',
        'lxml',
        'mp_api',
        'glom',
        'python-benedict[all]',
        # 'phonopy',
    ],
    entry_points={
        'console_scripts': [
            'fpflow=fpflow.scripts.fpflow:fpflow',
        ],
    },
)
#endregion

#region classes
#endregion