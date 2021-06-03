#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a standalone version of the example/index.html file by inlining all
javascript/css.

This creates a file that can be distributed in supplementary information
"""
import os
import base64
import argparse

from bs4 import BeautifulSoup
import requests


ROOT = os.path.join(os.path.dirname(__file__), "..")


def get_html():
    with open(os.path.join(ROOT, "app", "index.html")) as fd:
        html = fd.read()
    return BeautifulSoup(html, "html.parser")


def get_url(url):
    if url.startswith("http://") or url.startswith("https://"):
        data = requests.get(url)
        return data.text
    else:
        if url == "chemiscope-app.min.js":
            root = os.path.join(ROOT, "dist")
        else:
            root = os.path.join(ROOT, "app")

        path = os.path.realpath(os.path.join(root, url))
        if not os.path.exists(path):
            raise Exception(f"{path} does not exists, did you run `npm run build`?")

        with open(path) as fd:
            return fd.read()


def inline_ressources(html):
    for script in html.find_all("script"):
        if script.has_attr("src"):
            url = script["src"]
            if url.startswith("https://") or url.startswith("http://"):
                continue

            content = get_url(url)
            src = "data:text/javascript;base64,{}".format(
                base64.b64encode(content.encode("utf8")).decode("utf8")
            )

            new_script = html.new_tag("script", src=src)
            if script.has_attr("defer"):
                new_script.attrs["defer"] = None
            new_script.string = content
            script.replace_with(new_script)

    for link in html.find_all("link"):
        if link.attrs.get("rel", "") == ["stylesheet"]:
            url = link["href"]
            if url.startswith("https://") or url.startswith("http://"):
                continue
            content = get_url(url)
            style = html.new_tag("style")
            style.string = content
            link.replace_with(style)


def main():
    """
    Command-line utility to generate an stand-alone HTML file containing the
    chemiscope app. The standalone HTML can be loaded in a browser, and used to
    read chemiscope JSON inputs. In addition, a chemiscope JSON can be appended
    to the HTML to have a fully self-contained viewer for a chosen data set.
    """

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument(
        "input",
        type=str,
        default="",
        nargs="?",
        help="chemiscope JSON input that is appended to provide default data",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="output file name. Defaults to chemiscope.html or the name of the JSON",
    )
    args = parser.parse_args()

    html = get_html()
    inline_ressources(html)

    output_name = args.output
    if output_name == "":
        if args.input == "":
            output_name = "chemiscope.html"
        else:
            output_name = os.path.splitext(args.input)[0] + ".html"

    json_data = ""
    if args.input != "":
        if args.input.endswith(".gz"):
            raise ValueError(
                "Only plain-text JSON inputs can be used with the standalone viewer"
            )
        with open(args.input) as fd:
            json_data = fd.read()

    with open(output_name, "w") as fd:
        fd.write(str(html))
        fd.write(
            """<script type="text/javascript">
            const toHide = document.getElementsByClassName("hide-if-standalone");
            for (const element of toHide) {
                element.style.display = "none";
            }
            </script>"""
        )
        fd.write('<script id=standalone-json-data type="application/json">')
        fd.write(json_data)

    print("created " + output_name + "\nHappy chemiscoping!")


if __name__ == "__main__":
    main()
