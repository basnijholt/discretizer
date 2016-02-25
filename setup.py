from setuptools import setup

setup(name='discretizer',
      version='0.2',
      description='discretization module for kwant',
      url='https://gitlab.kwant-project.org/r-j-skolasinski/discretizer',
      author='Rafal Skolasinski and Sebastian Rubbert',
    #   author_email='',
      license='BSD 2-clause',
      packages=['discretizer'],
      install_requires=['sympy', 'numpy'],
      zip_safe=False)
