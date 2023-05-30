.. include:: links.inc

Release Tasks
=============

#. Setup `pypi`_ and `testpypi`_ accounts if you do not already have ones. Then
   put the following into your ``~/.pypirc``::

     [distutils]
         pypi
         testpypi

     [pypi]
     username = your_pypi_name
     password = your_pypi_password

     [testpypi]
     repository = https://test.pypi.org/legacy/
     username = your_testpypi_name
     password = your_testpypi_password

   Then register your accounts::

     # testpypi
     python setup.py register -r https://testpypi.python.org/pypi

     # pypi
     python setup.py register

   To be able to upload to pypi/testpypi without sending your password in plain
   text, install `twine`_. An upload is then done simply by::

     # testpypi
     twine upload -r testpypi <filenames>

     # pypi
     twine upload <filenames> # Be careful!

   **WARNING:** Uploads to `pypi`_ are permanent and cannot be updated later! If
   a wrong file is uploaded, the only fix is to create a new release. Always
   use `testpypi`_ first.

#. Bump version number:

   - in ``setup.py``
   - in ``docs/conf.py``

#. Regenerate and review the documentation::

     python setup.py build_ext -i build_sphinx
     firefox build/sphinx/html/index.html

#. Install and test the git version that is to be released. In the git
   repository, do::

     python setup.py install --user
     cd ..
     nosetests -v scikits.umfpack
     cd scikit-umfpack

#. Create a source distribution tarball, install it and test it::

     python setup.py sdist
     # unpack it, cd to it, then:
     python setup.py install --user
     cd ..
     nosetests -v scikits.umfpack

#. If OK, merge the version branch.

#. Upload to `testpypi`_ and test::

     python setup.py sdist
     twine upload -r testpypi dist/scikit-umfpack-<version>.tar.gz

   Note: if the upload fails with `This filename has previously been used, you
   should use a different version.`, just change the version in ``setup.py``
   to another value.

   Create a test install using `virtualenv`_::

     cd tmp/

     unset PYTHONPATH

     virtualenv test-umfpack

     cd test-umfpack
     source bin/activate

     python -m pip install --upgrade pip
     pip install nose numpy scipy
     pip install --no-deps --upgrade --pre -i https://testpypi.python.org/pypi scikit-umfpack

     nosetests -v scikits.umfpack

     deactivate

#. If OK, tag the version in git & push to github.

#. Upload to `pypi`_:

   - Check the version numbers in ``setup.py`` and ``docs/conf.py``
   - Change ISRELEASED in ``setup.py`` to True
   - Do::


       python setup.py sdist
       twine upload dist/scikit-umfpack-<version>.tar.gz

   - Change ISRELEASED in ``setup.py`` to False
   - For testing, see the previous step.

#. Update gh-pages::

     ./docs/do-gh-pages.sh
     git push -f origin gh-pages

#. If wheels are available for the released version, upload them also using
   `twine`_.

#. The wheels from https://github.com/scikit-umfpack/scikit-umfpack-wheels can
   be uploaded by
   https://github.com/MacPython/terryfy/blob/master/wheel-uploader script -
   save it to the scikit-umfpack directory, and change permissions to make it
   executable. Then, for example, upload to `testpypi`_::

     VERSION=0.3.1
     CDN_URL=https://3f23b170c54c2533c070-1c8a9b3114517dc5fe17b7c3f8c63a43.ssl.cf2.rackcdn.com
     ./wheel-uploader -u $CDN_URL -v -w ./wheels -t all -r testpypi scikit_umfpack $VERSION


   and to `pypi`_::

     ./wheel-uploader -u $CDN_URL -v -w ./wheels -t all scikit_umfpack $VERSION
