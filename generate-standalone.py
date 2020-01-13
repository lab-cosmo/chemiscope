#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate a standalone version of the example/index.html file by inlining all
javascript/css except for jsmol.

This creates a file that can be distributed in supplementary information
"""
import os

from bs4 import BeautifulSoup
import requests


ROOT = os.path.dirname(__file__)


def get_html():
    with open(os.path.join(ROOT, "examples", "standalone.in.html")) as fd:
        html = fd.read()
    return BeautifulSoup(html, 'html.parser')


def get_url(url):
    if url.startswith('http://') or url.startswith('https://'):
        data = requests.get(url)
        return data.text
    else:
        if url == 'sketchviz.min.js':
            path = os.path.join(ROOT, "dist")
        else:
            path = os.path.join(ROOT, "examples")

        with open(os.path.join(path, url)) as fd:
            return fd.read()


def inline_ressources(html):
    for script in html.find_all('script'):
        if script.has_attr('src'):
            url = script['src']
            if url.startswith('https://chemapps.stolaf.edu'):
                # don't replace jmol
                continue

            content = get_url(url)

            new_script = html.new_tag('script', type="text/javascript")
            new_script.string = content
            script.replace_with(new_script)


def main():
    html = get_html()
    inline_ressources(html)

    with open('standalone.html', 'w') as fd:
        fd.write(str(html))
        fd.write('<script id=json-data type="application/json">')

    print('created standalone.html!')


if __name__ == '__main__':
    main()
