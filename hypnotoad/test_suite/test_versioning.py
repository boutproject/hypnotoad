from hypnotoad.__version__ import get_versions


class TestVersioning:
    def test_init_version(self):
        from hypnotoad import __version__

        assert __version__ == get_versions()["version"]

    def test_git_hash(self):
        # Try to get git hash directly instead of through versioneer's get_versions()

        from boututils.run_wrapper import shell, shell_safe

        # check if git exists
        retval, git_version = shell("git --version")
        if retval == 0:
            # git exists

            from pathlib import Path
            from hypnotoad.__init__ import __file__ as hypnotoad_init_file

            hypnotoad_path = Path(hypnotoad_init_file).parent

            # check if hypnotoad is in it's own git repo. hypnotoad.__init__.py
            # should be in the hypnotoad/ subdirectory of the git repo if it is.
            # So the parent directory of hypnotoad_path should contain a '.git'
            # directory if hypnotoad is in a git repo
            if hypnotoad_path.parent.joinpath(".git").is_dir():
                retval, git_hash = shell_safe(
                    "cd "
                    + str(hypnotoad_path)
                    + '&& git describe --always --abbrev=0 --dirty --match "NOT A TAG"',
                    pipe=True,
                )
                git_hash = git_hash.strip()

                # found a git hash, check it is consistent with the one from versioneer
                dirty_pos = git_hash.find("-dirty")
                if dirty_pos != -1:
                    # Check versioneer says the repo is dirty
                    assert get_versions()["dirty"]

                    # remove "-dirty" from git_hash
                    git_hash = git_hash[:dirty_pos]
                else:
                    # repo is clean
                    assert not get_versions()["dirty"]

                assert git_hash == get_versions()["full-revisionid"]
        elif retval == 127:
            # git not installed
            pass
        else:
            raise RuntimeError(
                "'git --version' failed with {retval}. Output was {git_version}"
            )
