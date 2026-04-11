class Clann < Formula
  desc "Investigating phylogenetic information through supertree analyses"
  homepage "https://github.com/ChrisCreevey/clann"
  url "https://github.com/ChrisCreevey/clann/archive/refs/tags/V5.0.tar.gz"
  sha256 "a047e592b301aef8646f4284832ce8a8caa789e2927d9cf7678df04f2d8523eb"
  version "5.0.0"
  license "GPL-2.0-or-later"

  head "https://github.com/ChrisCreevey/clann.git", branch: "master"

  # GNU readline (the Homebrew version) gives full tab-completion and
  # history support.  Apple's libedit is used as a fallback by configure
  # but produces inferior interactive behaviour.
  depends_on "readline"

  # OpenMP parallel search replicates.
  # configure.ac auto-detects Apple Clang and switches to Homebrew GCC for
  # OpenMP support.  Declaring gcc here ensures it is installed so that
  # auto-detection succeeds on macOS.
  on_macos do
    depends_on "gcc"
  end

  def install
    system "./configure", "--disable-silent-rules",
                          "--prefix=#{prefix}",
                          "--with-readline=#{Formula["readline"].opt_prefix}"
    system "make", "-j#{ENV.make_jobs}"
    system "make", "install"
  end

  test do
    # A bare --help call should print usage and exit 0.
    assert_match "supertree", shell_output("#{bin}/clann --help 2>&1")

    # A minimal heuristic search on a tiny inline tree file should produce
    # a best tree line.
    (testpath/"test.ph").write <<~EOS
      (A,(B,(C,D)));
      (A,(C,(B,D)));
      (A,(D,(B,C)));
    EOS
    output = shell_output("#{bin}/clann hs #{testpath}/test.ph nreps=1 2>&1")
    assert_match "Supertree", output
  end
end
