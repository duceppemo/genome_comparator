# Contributing to genome_comparator

Thanks for your interest! Bug reports, questions, ideas and pull requests are all welcome.
Please follow the [code of conduct](CODE_OF_CONDUCT.md).

## Reporting a bug or asking a question
[Open an issue](https://github.com/duceppemo/genome_comparator/issues/new/choose) and pick the matching template.
For bugs, please include the version (`genome-comparator --version`), the command you ran and the log file
(`<output>/genome_comparator.log`, ideally from a run with `-v`). Check the
[troubleshooting page](https://github.com/duceppemo/genome_comparator/wiki/Troubleshooting) and the
[FAQ](https://github.com/duceppemo/genome_comparator/wiki/FAQ) first.

## Development setup
```
git clone https://github.com/duceppemo/genome_comparator
cd genome_comparator
conda env create -f environment.yml
conda activate genome_comparator
pip install --no-deps -e .
pytest --cov
```
The tests run on GitHub Actions for Python 3.10, 3.12 and 3.14, and coverage is reported on
[Codecov](https://app.codecov.io/gh/duceppemo/genome_comparator).

## Pull requests
* Branch from `master` and keep each pull request focused on one change.
* Add or update tests for any change in behaviour.
* Match the style of the surrounding code.
* Update the documentation in `wiki/` if the change affects users (see below), and add a line to
  `wiki/Changelog.md` under "Unreleased".

## Documentation
The [wiki](https://github.com/duceppemo/genome_comparator/wiki) is maintained in the [`wiki/`](wiki/) folder of this
repository and published automatically when changes reach `master`. **Do not edit the wiki on GitHub directly**:
those changes would be overwritten. Edit the files in `wiki/` in your pull request instead.

## Releases (maintainers)
1. Update the version in `genome_comparator/__init__.py` and `CITATION.cff` (`version` and `date-released`).
2. Rename "Unreleased" in `wiki/Changelog.md` to the version and date.
3. Commit, tag (`git tag -a vX.Y.Z`), push the commit and the tag, then create the GitHub release once the tests pass.
4. Publishing the GitHub release automatically:
   * uploads the package to [PyPI](https://pypi.org/project/genome-comparator/) (`.github/workflows/publish.yml`,
     trusted publishing: no token needed; the release tag must match the package version);
   * archives it on Zenodo with a new DOI (this can take a while when Zenodo is busy). Add the version DOI to the
     `identifiers` in `CITATION.cff`; the concept DOI in the README never changes.
5. Bioconda: update `recipe/meta.yaml` (version and sha256 of the new GitHub tarball). Bioconda's bot usually opens
   the pull request to `bioconda-recipes` on its own for new releases; otherwise copy the recipe and open it by hand.
