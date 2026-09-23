# Contributing

Bug reports, questions and pull requests are welcome. See
[CONTRIBUTING.md](https://github.com/duceppemo/genome_comparator/blob/master/CONTRIBUTING.md) for how to report a bug,
set up a development environment, open a pull request and make a release, and the
[code of conduct](https://github.com/duceppemo/genome_comparator/blob/master/CODE_OF_CONDUCT.md).

## Continuous integration
* Tests run on Python 3.10, 3.12 and 3.14 for every push and pull request.
* Coverage is uploaded to Codecov from the Python 3.12 job with the `CODECOV_TOKEN` secret. A failed upload fails
  the build on pushes only; pull requests from forks upload without the token, and Dependabot pull requests (which
  cannot read secrets) skip the upload.
* Dependabot opens one pull request per month to update the GitHub Actions used by the workflows
  (`.github/dependabot.yml`).

## Editing this wiki
This wiki is maintained in the [`wiki/`](https://github.com/duceppemo/genome_comparator/tree/master/wiki) folder of
the main repository and published to the GitHub wiki automatically when changes are pushed to `master`.
**Do not edit the wiki on GitHub directly**: changes would be overwritten. Edit the files in `wiki/` and open a pull
request instead.

* Page file names become page titles: `Output-files.md` → "Output files".
* Link to other pages without the extension: `[Usage](Usage)`.
* `_Sidebar.md` is the navigation shown on every page.
* Images are stored in the repository's `assets/` folder and linked with their
  `https://raw.githubusercontent.com/duceppemo/genome_comparator/master/assets/...` URL.
