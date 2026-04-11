# Packaging Clann

This document describes how to install Clann via Conda or Homebrew, and how to
submit Clann to the Bioconda and Homebrew-bio community channels.

---

## Conda (Bioconda)

### Installing from Bioconda (once the recipe is accepted)

```bash
conda install -c bioconda -c conda-forge clann
```

Or with `mamba` for faster solves:

```bash
mamba install -c bioconda -c conda-forge clann
```

### Building locally from the recipe

Requirements: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or
Miniforge, with `conda-build` installed.

```bash
# Install build tools (once)
conda install -y conda-build conda-verify

# Build the package (from the clann repository root)
conda build conda-recipe/

# Install the locally built package
conda install --use-local clann
```

The recipe lives in `conda-recipe/` and follows the
[Bioconda packaging guidelines](https://bioconda.github.io/contributor/guidelines.html).

### Submitting to Bioconda

1. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
   on GitHub.
2. Copy `conda-recipe/` into `recipes/clann/` in your fork:

   ```bash
   cp -r conda-recipe/ path/to/bioconda-recipes/recipes/clann/
   ```

3. Lint the recipe:

   ```bash
   cd path/to/bioconda-recipes
   bioconda-utils lint --packages clann
   ```

4. Open a Pull Request against `bioconda/bioconda-recipes`. The Bioconda CI
   bot will build and test the recipe on Linux and macOS.

### Updating the recipe for a new release

When a new tag (e.g. `V5.1`) is pushed to GitHub:

1. Update `{% set version %}` and `{% set tag %}` in `conda-recipe/meta.yaml`.
2. Update `sha256` with the checksum of the new tarball:

   ```bash
   curl -sL https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz \
     | sha256sum
   ```

3. Reset `build: number` to `0`.
4. Open a PR to bioconda-recipes with the updated recipe.

---

## Homebrew

### Installing from a tap (once published)

```bash
# Add the tap (replace <user> with the GitHub account hosting the tap)
brew tap ChrisCreevey/clann
brew install clann
```

Or, once accepted into [homebrew-bio](https://github.com/brewsci/homebrew-bio):

```bash
brew tap brewsci/bio
brew install clann
```

### Using the formula directly (without a tap)

```bash
brew install --formula Formula/clann.rb
```

Or install straight from the repository URL:

```bash
brew install https://raw.githubusercontent.com/ChrisCreevey/clann/master/Formula/clann.rb
```

### Setting up a personal Homebrew tap

1. Create a new GitHub repository named **`homebrew-clann`** under your
   account (the `homebrew-` prefix is required by Homebrew).
2. Create a `Formula/` directory inside it and copy `Formula/clann.rb` there.
3. Users can then install with:

   ```bash
   brew tap ChrisCreevey/clann
   brew install clann
   ```

### Submitting to homebrew-bio

[brewsci/homebrew-bio](https://github.com/brewsci/homebrew-bio) is the
community Homebrew tap for bioinformatics tools, analogous to Bioconda.

1. Fork `brewsci/homebrew-bio`.
2. Copy `Formula/clann.rb` into the `Formula/` directory of your fork.
3. Verify it builds cleanly:

   ```bash
   brew install --build-from-source Formula/clann.rb
   brew test Formula/clann.rb
   brew audit --strict Formula/clann.rb
   ```

4. Open a Pull Request against `brewsci/homebrew-bio`.

### Updating the formula for a new release

1. Update `url` to the new tag (e.g. `V5.1`).
2. Update `sha256`:

   ```bash
   curl -sL https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz \
     | sha256sum
   ```

3. Update `version`.
4. Open a PR to the relevant tap.

---

## GitHub Releases

Both Bioconda and Homebrew reference the source tarball that GitHub generates
automatically from a tag.  No manual upload is needed — pushing a tag is
sufficient:

```bash
git tag V5.1
git push origin V5.1
```

GitHub will serve the tarball at:

```
https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz
```

---

## Continuous Integration

A GitHub Actions workflow (`.github/workflows/conda-build.yml`) builds and
tests the Conda package on Ubuntu and macOS whenever `conda-recipe/` is
changed or a new version tag is pushed.  This gives confidence that the recipe
is correct before submitting to Bioconda.
