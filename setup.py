
import os

from setuptools import setup

setup (
    name='Seminario',
    version='1.0.0',
    description='Tool to use the Seminario method',
    long_description=open('README.md').read(),
    url='https://github.com/sgenheden/Seminario',
    author='Samuel Genheden',
    author_email='samuel.genheden@gmail.com',
    license='GNU General Public Licence',

    packages=['seminario',],
    entry_points={'console_scripts': ['seminario_ff = seminario.tools:seminario_ff', ]},
    install_requires=['parmed','numpy',],
)
