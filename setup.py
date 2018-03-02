import setuptools
import sys

if sys.version_info < (3, 0):
    raise EnvironmentError('Please install using pip3 or python3')

setuptools.setup(author='Chris Rosenthal',
                 author_email='crosenth@gmail.com',
                 description='Ribosomal RNA Allele Tree',
                 keywords=['ncbi', 'rrndb', 'genetics', 'genomics'],
                 name='rrat',
                 packages=setuptools.find_packages(exclude=['tests']),
                 entry_points={
                     'console_scripts': {'rrat = rrat:main'}},
                 version=0.1,
                 url='https://github.com/crosenth/rrat',
                 package_data={'classifier': ['data/*']},
                 license='GPLv3',
                 classifiers=[
                     'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                     'Development Status :: 4 - Beta',
                     'Environment :: Console',
                     'Operating System :: OS Independent',
                     'Intended Audience :: End Users/Desktop',
                     'License :: OSI Approved :: '
                     'GNU General Public License v3 (GPLv3)',
                     'Programming Language :: Python :: 3 :: Only'])
