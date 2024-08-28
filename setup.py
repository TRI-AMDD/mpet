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
    url='https://github.com/TRI-AMDD/mpet',
    packages=[
        'mpet','mpet.plot',
        'mpet.electrode.diffusion',
        'mpet.electrode.materials',
        'mpet.electrode.reactions',
        'mpet.electrolyte','mpet.config'
    ],
    install_requires=['numpy','scipy','matplotlib','pyQt5', 'h5py', 'configparser', 'schema',
                      'daetools @ https://sourceforge.net/projects/daetools/files/daetools/2.3.0/daetools-2.3.0-gnu_linux-x86_64.zip ; python_version >= "3.10" and python_version <= "3.12" and sys_platform == "linux"',
                      'daetools @ https://sourceforge.net/projects/daetools/files/daetools/2.3.0/daetools-2.3.0-win64.zip ; python_version >= "3.10" and python_version <= "3.12" and sys_platform == "win32"',
                      'daetools @ https://sourceforge.net/projects/daetools/files/Old-releases/daetools-old/1.9.0/daetools-1.9.0-gnu_linux-x86_64.tar.gz ; python_version >= "3.5" and python_version <= "3.7" and sys_platform == "linux"',
                      'daetools @ https://sourceforge.net/projects/daetools/files/Old-releases/daetools-old/1.9.0/daetools-1.9.0-win64.zip ; python_version >= "3.5" and python_version <= "3.7" and sys_platform == "win32"',
                      ],
    extras_require={'test':['pytest','coverage', 'coveralls', 'flake8'],
                    'doc':['sphinx','sphinx_rtd_theme'],
                    'dashboard': ['dash', 'dash_bootstrap_components'],
                    'cluster_jobs': ['dask-jobqueue', 'bokeh']},
    python_requires='>=3.6',
    scripts=['bin/mpetrun.py','bin/mpetplot.py','bin/run_jobs.py', 'bin/create_ensemble.py',
             'bin/mpet_create_runjobs_dashboard.py', 'bin/mpet_plot_app.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
