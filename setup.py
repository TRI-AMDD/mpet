import setuptools

exec(open('mpet/version.py').read())

with open("README.md", "r", encoding="utf-8") as fh:
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
    packages=['mpet','mpet.plot','mpet.electrode'],
    install_requires=["numpy","scipy","matplotlib","pyQt5"],
    extras_require = {'test':['pytest','coverage', 'coveralls','configparser','h5py']},
    python_requires='>=3.5,<3.8',
    scripts=['bin/mpetrun.py','bin/mpetplot.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)
