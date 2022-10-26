import setuptools
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

here = os.path.abspath(os.path.dirname(__file__))

about = {}
with open(os.path.join(here, 'src', '__version__.py'), mode='r',
          encoding='utf-8') as f:
    exec(f.read(), about)

with open('requirements.txt') as fp:
    install_requires = fp.read()

setuptools.setup(
    name=about['__title__'],
    version=about['__version__'],
    author="XYZ",
    author_email="xyz@abc.com",
    description="Open-source DBS",
    packages=setuptools.find_packages(),
    python_requires='>=3.8',
    install_requires=install_requires,
)
