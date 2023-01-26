Release process
---------------

1. Pick the new version number, e.g. 0.1.2, for the next release (bumping
   major, minor, or patch version according to semantic versioning)
2. Ensure doc/whats-new.md is updated with the number of the version to be
   released.
3. Create a new tag and release with the 'Create a new release' link on the
   home page. Name the new release according to 2., e.g. 0.1.2 with no prefix.
4. Copy the section for this version from doc/whats-new.md into the release
   notes.
5. Publish the release! The release will automatically upload hypnotoad to
   PyPi.
6. Approve the automatic PR from conda-forge to release the new version on
   conda.
7. Automatic Zenodo integration with the Github repo has been activated. For
   the next version after 0.5.1, the Zenodo entry should be automatically
   created. Check that the automatic version links to the same record as the
   [manually created ones](https://doi.org/10.5281/zenodo.6360326) and that the
   metadata updated during the manual process was updated correctly.

   The 'manual' process was:

   Make a new Zenodo entry for the new version on
   https://doi.org/10.5281/zenodo.6360326.  

   i. Upload the .tar.gz from the GitHub release.  
   ii. Update the 'publication date'.  
   iii. Check for any new contributors since the previous release.  
   iv. Update the 'Version' tag.  
   v. Update the link in 'Related/alternate identifiers' - the version number
      will need changing.  
   vi. Save and publish the release!
