import setuptools

version = {}
with open('mpet/version.py', 'r', encoding='utf-8') as fh:
    exec(fh.read(), version)
__version__ = version['__version__']

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='mpet',
    version=__version__,
    description='Multiphase porous electrode theory',
    long_description=long_description,
    long_description_content_type="text/markdown",
    author='Dan Cogswell',
    author_email='cogswell@mit.edu',
    license='MIT',
    url='https://bitbucket.org/bazantgroup/mpet',
    packages=['mpet','mpet.plot','mpet.electrode','mpet.config'],
    install_requires=['numpy','scipy','matplotlib','pyQt5', 'h5py', 'configparser', 'schema', 'dask-jobqueue'],
    extras_require={'test':['pytest','coverage', 'coveralls', 'flake8'],
                    'doc':['sphinx','sphinx_rtd_theme']},
    python_requires='>=3.5',
    scripts=['bin/mpetrun.py','bin/mpetplot.py','bin/run_jobs.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
