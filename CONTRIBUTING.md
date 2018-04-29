Contributing to RichDEM
========================

Project Description
-------------------

RichDEM is a set of digital elevation model (DEM) hydrologic analysis tools.
RichDEM uses parallel processing and state of the art algorithms to quickly
process even very large DEMs.

RichDEM offers a variety of flow metrics, such as D8 and D∞. It can flood or
breach depressions. It can calculate flow accumulation, slops, curvatures, &c.

RichDEM is available as a performant C++ library, a low-dependency Python
package, and a set of command-line tools.

- https://richdem.readthedocs.io/en/latest/



Contributor Agreement
---------------------

By contributing to this project you agree to the following Developer Certificate
of Origin:

    Developer Certificate of Origin
    Version 1.1

    Copyright (C) 2004, 2006 The Linux Foundation and its contributors.
    1 Letterman Drive
    Suite D4700
    San Francisco, CA, 94129

    Everyone is permitted to copy and distribute verbatim copies of this
    license document, but changing it is not allowed.


    Developer's Certificate of Origin 1.1

    By making a contribution to this project, I certify that:

    (a) The contribution was created in whole or in part by me and I
        have the right to submit it under the open source license
        indicated in the file; or

    (b) The contribution is based upon previous work that, to the best
        of my knowledge, is covered under an appropriate open source
        license and I have the right under that license to submit that
        work with modifications, whether created in whole or in part
        by me, under the same open source license (unless I am
        permitted to submit under a different license), as indicated
        in the file; or

    (c) The contribution was provided directly to me by some other
        person who certified (a), (b) or (c) and I have not modified
        it.

    (d) I understand and agree that this project and the contribution
        are public and that a record of the contribution (including all
        personal information I submit with it, including my sign-off) is
        maintained indefinitely and may be redistributed consistent with
        this project or the open source license(s) involved.



Developer Resources
-------------------

RichDEM code is hosted on GitHub:

* https://github.com/r-barnes/richdem



Issue Tracking
--------------

RichDEM uses Github to track ongoing development and issues:

* https://github.com/r-barnes/richdem/issues



Building
--------

RichDEM requires a C++ compiler to build. It is regularly tested with GCC and
clang. It is not regularly tested with MS Visual Studio.



Contributing
------------

RichDEM uses git pull requests for contributions. To create a pull request, follow these steps:

* Fork the RichDEM project on GitHub - go to https://github.com/r-barnes/richdem and click 'Fork'.
* Create a branch on your forked project that contains your work. See 'Coding Standards', below.
* Use GitHub to open a pull request against the r-barnes RichDEM repository - from your branch on
  GitHub, click 'New Pull Request'.
* Respond to comments on your pull request as they are made.
* When ready, your pull request will be merged by an official RichDEM contributor.



Coding Standards
----------------

* An initial pull request should consist of a single commit, rebased on top of
  the current master branch.
  * Additional commits can be added in response to code review comments. Github 
    automagically updates the pull request if these commits are made to the
    feature branch mentioned above.
* The commit message must be descriptive.
* Code should conform to the style guide below.
* Code should include unit tests when appropriate.



License and Copyright
---------------------

RichDEM is provided under the GPL license and any contributions must maintain
this. 

In the future: To ensure proper licensing, source files must contain an
appropriate license header.



Contact
-------

Contact about bugs should be made via Github issues. If an issue is not acted on
within a reasonable timeframe, please email Richard.

Contact about collaborations, larger feature requests, and such can be made via
email to Richard. Please be mindful that a Github issue is preferable to an
email.
