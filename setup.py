from setuptools import setup, find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()

here = os.path.abspath(os.path.dirname(__file__))

about = {}
with open(os.path.join(here, 'ossdbs', '__version__.py'),
          mode='r',
          encoding='utf-8') as f:
    exec(f.read(), about)

with open('requirements.txt') as fp:
    install_requires = fp.read()

setup(
    name=about['__title__'],
    version=about['__version__'],
    author="XYZ",
    author_email="xyz@abc.com",
    description="Open-source DBS",
    packages=find_packages(include=['ossdbs', 'ossdbs.*']),
    python_requires='>=3.8',
    install_requires=install_requires,
    entry_points={'console_scripts': ['ossdbs=ossdbs.main:main']}
)
