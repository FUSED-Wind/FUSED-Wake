from setuptools import setup, find_packages  # Always prefer setuptools over distutils

setup(
    name='py4we',
    version='0.0.1a1',
    packages=find_packages(),
    url='https://github.com/piredtu/py4we',
    license='Apache 2.0',
    author='Pierre-Elouan Rethore',
    author_email='pire@dtu.dk',
    description='A pure python set of wrappers for wind energy tools',
    classifiers=[
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'Programming Language :: Python :: 2.7'],
)
