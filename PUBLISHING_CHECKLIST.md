# Maintainer Publishing Checklist

This document is a step-by-step guide for **Chris Creevey** to publish Clann
to the Homebrew and Bioconda packaging channels and to keep those packages
up to date on each new release.

The packaging files already live in this repository:

| File | Purpose |
|---|---|
| `Formula/clann.rb` | Homebrew formula |
| `conda-recipe/meta.yaml` | Bioconda recipe metadata |
| `conda-recipe/build.sh` | Bioconda build script |
| `.github/workflows/conda-build.yml` | CI that validates the conda recipe |

---

## One-time prerequisites

### Homebrew

- [ ] Install Homebrew if it is not already present:
  ```bash
  /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  ```
- [ ] Verify `brew` is in your `PATH`:
  ```bash
  brew --version
  ```

### Conda / Bioconda

- [ ] Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or
  [Miniforge](https://github.com/conda-forge/miniforge) if it is not already
  present.
- [ ] Install the build tools:
  ```bash
  conda install -y conda-build conda-verify
  ```
- [ ] Have a GitHub account with SSH or HTTPS access configured (needed for
  forking bioconda-recipes).

---

## Part 1 — Homebrew

There are two levels here: a **personal tap** (you control it, users install
from it today) and the **homebrew-bio community tap** (wider reach but
requires a PR review).  Do them in order.

### Step 1 — Create the `homebrew-clann` tap repository (first release only)

1. Go to <https://github.com/new>.
2. Set the repository name to exactly **`homebrew-clann`**.
   (The `homebrew-` prefix is required by Homebrew's tap convention.)
3. Leave it public, with no README.
4. Create it.
5. Clone it locally:
   ```bash
   git clone git@github.com:ChrisCreevey/homebrew-clann.git
   cd homebrew-clann
   ```
6. Create a `Formula/` directory and copy the formula into it:
   ```bash
   mkdir Formula
   cp /path/to/clann/Formula/clann.rb Formula/
   git add Formula/clann.rb
   git commit -m "Add clann formula"
   git push
   ```
7. Test that Homebrew can install from your tap:
   ```bash
   brew tap ChrisCreevey/clann
   brew install clann
   clann --help
   ```

From this point on, anyone can install Clann with:
```bash
brew tap ChrisCreevey/clann
brew install clann
```

---

### Step 2 — Submit to homebrew-bio (community tap)

[brewsci/homebrew-bio](https://github.com/brewsci/homebrew-bio) is the
community tap for bioinformatics tools.  Getting accepted here means users
can install via `brew tap brewsci/bio && brew install clann` without needing
your personal tap.

1. Fork `brewsci/homebrew-bio` on GitHub.
2. Clone your fork:
   ```bash
   git clone git@github.com:ChrisCreevey/homebrew-bio.git
   cd homebrew-bio
   ```
3. Copy the formula:
   ```bash
   cp /path/to/clann/Formula/clann.rb Formula/
   git add Formula/clann.rb
   git commit -m "Add clann formula"
   git push
   ```
4. Audit and test locally before opening the PR:
   ```bash
   brew install --build-from-source Formula/clann.rb
   brew test Formula/clann.rb
   brew audit --strict Formula/clann.rb
   ```
   All three commands must complete without errors.
5. Open a Pull Request from your fork to `brewsci/homebrew-bio`.
   - Title: `Add clann`
   - Description: brief summary (one sentence) and a link to the Clann GitHub page.
6. Wait for the CI checks on the PR to go green, then respond to any reviewer
   feedback.

---

## Part 2 — Bioconda

### Step 3 — Validate the recipe locally before submitting

```bash
cd /path/to/clann
conda build conda-recipe/
conda build --test conda-recipe/
```

Both commands must complete without errors on Linux (and ideally macOS).
The GitHub Actions workflow (`.github/workflows/conda-build.yml`) does this
automatically on every push to `conda-recipe/` — check that it is green before
proceeding.

### Step 4 — Submit the recipe to bioconda-recipes (first release only)

1. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
   on GitHub.
2. Clone your fork and create a feature branch:
   ```bash
   git clone git@github.com:ChrisCreevey/bioconda-recipes.git
   cd bioconda-recipes
   git checkout -b add-clann
   ```
3. Copy the recipe directory:
   ```bash
   cp -r /path/to/clann/conda-recipe/ recipes/clann/
   git add recipes/clann/
   git commit -m "Add clann recipe"
   git push -u origin add-clann
   ```
4. Lint the recipe (install `bioconda-utils` first if needed):
   ```bash
   bioconda-utils lint --packages clann
   ```
   Fix any warnings or errors before opening the PR.
5. Open a Pull Request from your `add-clann` branch to `bioconda/bioconda-recipes`.
   - Title: `Add clann`
   - The Bioconda CI bot will build and test on Linux and macOS automatically.
6. Respond to reviewer feedback and wait for the PR to be merged.

Once merged, users can install with:
```bash
conda install -c bioconda -c conda-forge clann
```

---

## Part 3 — Releasing a new version

Do these steps every time you cut a new release (e.g. `V5.1`).

### Step 5 — Tag and push to GitHub

```bash
git tag V5.1
git push origin V5.1
```

GitHub automatically creates the source tarball at:
```
https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz
```

### Step 6 — Get the sha256 checksum of the new tarball

```bash
# Linux:
curl -sL https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz \
  | sha256sum

# macOS:
curl -sL https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz \
  | shasum -a 256
```

Note the 64-character hex string — you will need it in the next two steps.

### Step 7 — Update the Homebrew formula

Edit `Formula/clann.rb` in **this** repository and in your `homebrew-clann`
tap repository:

```ruby
url "https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.1.tar.gz"
sha256 "<new sha256 from step 6>"
version "5.1.0"
```

Then push to your personal tap:
```bash
cd /path/to/homebrew-clann
# edit Formula/clann.rb
git add Formula/clann.rb
git commit -m "clann 5.1.0"
git push
```

Open a PR to `brewsci/homebrew-bio` with the same change (if you were
accepted there).

### Step 8 — Update the Bioconda recipe

Edit `conda-recipe/meta.yaml` in **this** repository:

```yaml
{% set version = "5.1.0" %}
{% set tag = "V5.1" %}
...
source:
  sha256: <new sha256 from step 6>
build:
  number: 0       # reset to 0 for every new version
```

Commit and push the update:
```bash
git add conda-recipe/meta.yaml
git commit -m "Bump clann to 5.1.0"
git push
```

The GitHub Actions CI workflow will rebuild and test the package automatically.

Then open a PR to `bioconda/bioconda-recipes` with the same change to
`recipes/clann/meta.yaml`.

---

## Quick reference — files to edit per release

| File | Fields to update |
|---|---|
| `Formula/clann.rb` | `url` (tag), `sha256`, `version` |
| `conda-recipe/meta.yaml` | `version`, `tag`, `sha256`, reset `build: number` to `0` |
| `homebrew-clann/Formula/clann.rb` | same as `Formula/clann.rb` above |

---

## Useful links

| Resource | URL |
|---|---|
| bioconda-recipes | <https://github.com/bioconda/bioconda-recipes> |
| Bioconda contributor guide | <https://bioconda.github.io/contributor/guidelines.html> |
| homebrew-bio tap | <https://github.com/brewsci/homebrew-bio> |
| Homebrew formula cookbook | <https://docs.brew.sh/Formula-Cookbook> |
| Clann GitHub | <https://github.com/ChrisCreevey/clann> |
