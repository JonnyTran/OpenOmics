# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions
- Web development with the Dash for the OpenOmics dashboard webserver.
- Documentation standards for the OpenOmics Python API.
- Organize the library of genomics, functional ontologies, interactions, and sequence databases for variety of
  biological studies.
- Implement general purpose utilities for importing various fasta, gtf and sequencing files.

## Report Bugs

Report bugs at [openomics/issues](https://github.com/JonnyTran/openomics/issues).

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

## Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

## Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

## Write Documentation

OpenOmics could always use more documentation, whether as part of the
[official OpenOmics docs](https://openomics.readthedocs.io/), in docstrings within the API, or even on the web in blog posts,
articles, and such.

If you'd like to help write RTD documentations, note:
- Documentation pages are written in markdown using [myst-parser](https://myst-parser.readthedocs.io/en/latest/index.html)
- The Sphinx theme used is [furo](https://pradyunsg.me/furo/)
- The autodoc package used is [sphinx-automodapi](https://sphinx-automodapi.readthedocs.io/en/latest/)

## Submit Feedback

The best way to send feedback is to file an issue at [openomics/issues](https://github.com/JonnyTran/openomics/issues).

If you are proposing a feature:

- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.
- Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Get Started!

Ready to contribute? Here's how to set up `openomics` for local development.

1. Fork the `openomics` repo on GitHub.
2. Clone your fork locally and work on the develop branch:
```
$ git clone git@github.com:your_name_here/openomics.git
$ git checkout develop
```

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development:
```
$ mkvirtualenv openomics
$ cd openomics/
$ python setup.py develop
```

4. Create a branch for local development:
```
$ git checkout -b name-of-your-bugfix-or-feature
```
   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox:
```
$ flake8 openomics tests
$ python setup.py test or py.test $ tox
```

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub:
```
$ git add .
$ git commit -m "Your detailed description of your changes."
$ git push develop name-of-your-bugfix-or-feature
```
7. Submit a pull request through the GitHub website to the develop branch. Once major features are tested, we can create
   another pull-request to the **master** branch.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests. Run tests by with `pytest ./` and make sure tests are 100% passing.
2. If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a
   docstring, and add the feature to the list
   in [docs/history.md](https://github.com/JonnyTran/OpenOmics/blob/master/docs/history.md).
3. The pull request should work for Python 3.6 or higher, and for PyPi. Check
   [Github Actions Tests](https://github.com/JonnyTran/OpenOmics/actions/workflows/python-package.yml)
   and make sure that the tests pass for all supported Python versions and operating systems.

## Tips

To run the automated tests locally, run this at the root directory:

    pytest ./

```{hint}
To run a subset of tests:

    $ py.test tests.test_openomics

```

To run tests targeting various operating systems and Python versions, make a pull-request to the **master** branch which
will run as [Github Actions Tests](https://github.com/JonnyTran/OpenOmics/actions/workflows/python-package.yml).

## Deploying

A reminder for the maintainers on how to deploy. Make sure all your changes are committed (including an entry in
HISTORY.rst). Then run:

    $ bumpversion patch # possible: major / minor / patch
    $ git push --tags

Github Actions will then deploy to PyPI if tests pass.

## Code of Conduct
Please note that the OpenOmics project is released with a Contributor Code of Conduct. By contributing to this project you agree to abide by its terms.

[openomics/issues]: https://github.com/JonnyTran/openomics/issues
