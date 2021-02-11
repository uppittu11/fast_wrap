from setuptools import setup, find_packages
import sys

try:
    import mdtraj
except ImportError:
    print('Building and running analysis requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!')
    sys.exit(1)

setup(name='fast_wrap',
      version='0.1',
      description=('A fast script to wrap atoms or molecules into '
                   'a simulation box using MDTraj as a backend.'),
      url='https://github.com/uppittu11/fast_wrap',
      author='Parashara Shamaprasad',
      author_email='p.shama@vanderbilt.edu',
      license='MIT',
      packages=find_packages(),
      package_dir={'fast_wrap': 'fast_wrap'},
      include_package_data=False,
      install_requires=["mdtraj", "numpy"],
      #entry_points={
      #    "console_scripts" : ["analyze=bin.analyze:main"],
      #    },
      #scripts=['bin/progress.sh'],
)
