# Contributing guide

This guide will outline the basic process for contributing to the SNF-simulations project.

## Project management

Issues and pull requests are managed on [GitHub](https://github.com/ekneale/SNF-simulations).

**Before starting work:**
- Check the [GitHub Issues page](https://github.com/ekneale/SNF-simulations/issues) to see if the bug, feature, or enhancement you're interested in is already being discussed.
- If you can't find anything relevant, [raise a new issue](https://github.com/ekneale/SNF-simulations/issues/new) describing the change or problem. This helps coordinate effort and allows maintainers to provide early feedback.

**When opening a pull request:**
- Link your pull request to the relevant issue so context and discussion are kept together.
- Reference the issue number in your PR description (e.g., "Closes #123" or "Addresses #45") to auto-link the PR to the issue.

## Typical contribution workflow

Details for the development setup and tools are covered in the [development guide](development.md), but the general workflow for contributing is as follows:

1. Find or create an issue on GitHub as above, to track your work and any discussions.
2. Install the project in a local development environment.
3. Create a new branch, make and commit your changes.
4. Run local checks (e.g. `ruff`, `ty`, `pytest`).
5. Push your branch to GitHub and open a pull request, linking to the relevant issue.
6. Address any feedback and ensure all GitHub Actions checks pass.
7. Once approved, your PR will be merged and your contribution will be part of the next release.

## What to include in your branch

Here are some guidelines to ensure your contribution is complete and maintainable.

### Code

- Follow the existing code style and conventions used in the project.
- Include type annotations for any new functions or classes, following the existing style.
- Ensure your code is well-structured and modular, making it easier to test and maintain.

### Documentation
- Include docstrings for any new functions or classes, following the [Google style docstring format](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings).
  - This should ensure the new code is properly included in the [auto-generated API pages](../apidocs/index.rst).
- If your changes affect user-facing behavior or add new features, update the relevant documentation pages in [the user guide](../user_guide/index.md).
- If your changes affect the project installation, usage, or other high-level aspects, remember to update the package `README.md` file.

### Tests
- Add or update tests in `tests/` to cover your new code and any modified functionality.
- Run `pytest` with code coverage checks to verify that new additions are sufficiently tested.

## Licensing

This project is distributed under the BSD 3-Clause license. Please ensure contributions are compatible with this license.

When adding third-party code or data, include the appropriate attribution and license information in the repository documentation as needed. See `LICENSE.txt` for the full project license text.
