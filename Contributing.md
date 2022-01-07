# Contributing to chemfiles

:tada: First off, thanks for taking the time to contribute to chemiscope! :tada:

If you want to contribute but feel a bit lost, do not hesitate to contact us and
ask your questions! We will happily mentor you through your first contributions.

## Area of contributions

The first and best way to contribute to chemiscope is to use it and advertise it
to other potential users. Other than that, you can help with:

-   documentation: correcting typos, making various documentation clearer;
-   bug fixes and improvements to existing code;
-   and many more …

All these contributions are very welcome. We accept contributions via Github
pull request (have a look [here][pr] for Github model of pull request). If you
want to work on the code and pick something easy to get started, have a look at
the [first good issues][easy-issues].

## Bug reports and feature requests

Bug and feature requests should be reported as [Github issue][issue]. For bugs,
you should provide information so that we can reproduce it: what did you try?
What did you expect? What happened instead? Please provide any useful code
snippet or input file with your bug report.

If you want to add a new feature to chemiscope, please create an [issue] so that
we can discuss it, and you have more chances to see your changes incorporated.

### Code contribution check-list

Every item in this list is explained in the next section

-   [ ] Fork chemiscope;
-   [ ] Create a local branch;
-   [ ] Add code / correct typos / ...;
-   [ ] Check that the code passes lint checks;
-   [ ] Push to Github;
-   [ ] Create a Pull Request;
-   [ ] Discuss your changes with the reviewers;
-   [ ] Have your code merged
-   [ ] Celebrate! :tada: :cake: :tada:

### Contribution tutorial

In this small tutorial, you should replace `<angle brackets>` as needed. If
anything is unclear, please ask for clarifications! There are no dumb questions.

---

Start by [forking chemiscope][fork], and then clone and start a development
server for your fork by running:

```bash
git clone https://github.com/<YOUR USERNAME>/chemiscope
cd chemiscope
npm install
npm start
```

Then create a new branch for your changes

```bash
git checkout -b <new-branch>
```

Implement your changes, including the documentation (in `docs`) for new code.
You can then navigate to `http://localhost:8080` to interact with the code and
check that everything works as expected. It is always a good idea to check that
your code is working with multiple browsers (Chrome, Firefox, Edge, Safari).

Run tests and lints (we use [eslint] and [prettier] to ensure a consistent
coding style):

```bash
npm test
```

We suggest that you configure your code editor to automatically re-format the
code when you save the files. There are prettier plugins for most editors, see
the "Editor Support" section of the [prettier] website. If you want to manually
format your files, you can use

```bash
npx prettier --write <path/to/the/files>
```

Finally, you can push your code to Github, and create a [Pull Request][pr] to
the `lab-cosmo/chemiscope` repository.

```bash
git commit  # ask for help if you don't know how to use git
git push -u origin <new-branch>
```

[pr]: https://help.github.com/articles/using-pull-requests/
[easy-issues]: https://github.com/lab-cosmo/chemiscope/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22
[fork]: https://help.github.com/articles/fork-a-repo/
[issue]: https://github.com/lab-cosmo/chemiscope/issues/new
[eslint]: https://eslint.org/
[prettier]: https://prettier.io/
