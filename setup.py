import setuptools

exec(open('mpet/version.py').read())

setuptools.setup(
    name='mpet',
    version=__version__,
    description='Multiphase porous electrode theory',
    author='Dan Cogswell',
    author_email='cogswell@mit.edu',
    license='MIT',
    url='https://bitbucket.org/bazantgroup/mpet',
    packages=['mpet','mpet.plot','mpet.electrode'],
    install_requires=["numpy","scipy","matplotlib","pyQt5"],
    python_requires='>=3.5,<3.8',
    scripts=['bin/mpetrun.py','bin/mpetplot.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
)