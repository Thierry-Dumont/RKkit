from setuptools import setup
setup(
    name='RKkit',
    version='0.1.0',
    description='Playing with Runge-Kutta methods',
    url='https://github.com/Thierry-Dumont/RKkit',
    author='Thierry Dumont',
    author_email='tdumont@math.cnrs.fr',
    license='GNU GENERAL PUBLIC LICENSE Version 3',
    packages=['RKkit'],
    install_requires=[],
    extras_require={"passagemath":["passagemath-standard","passagemath-plot",
                                   "passagemath-symbolics",
                                   "passagemath-categories",
                                   "passagemath-combinat",
                                   "passagemath-flint",
                                   "passagemath-environment",
                                   "passagemath-singular",],},
)
