import setuptools

exec(open('mpet/version.py').read())

setuptools.setup(
    name='mpet',
    version=__version__,
    description='Multiphase porous electrode theory',
    author='Dan Cogswell',
    author_email='cogswell@mit.edu',
    url='https://bitbucket.org/bazantgroup/mpet',
    packages=['mpet'],
    python_requires='>=3.5,<3.7',
    scripts=['bin/mpetrun.py','bin/mpetplot.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)