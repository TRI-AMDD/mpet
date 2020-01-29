import setuptools

#exec(open('UI/version.py').read())

setuptools.setup(
    name='mpetUI',
    version='1.0',
    description='User Interface for MPET using VKML',
    author='Surya Mitra, Edwin García, Alex Bartol, Luke Robinson, Jarrod Lund',
    author_email='sayalaso@purdue.edu,redwing@purdue.edu',
    packages=['VKML','VKML.gui','VKML.gui.tk','UI'],
    scripts=['bin/mpetUI','bin/mpetviz'],
    classifiers=[
        "Programming Language :: Python :: 3.6",
    ],
    data_files=[('man/man1',['mpetUI.1'])]
)
