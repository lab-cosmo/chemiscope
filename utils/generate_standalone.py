#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Generate a standalone version of the example/index.html file by inlining all
javascript/css except for jsmol.

This creates a file that can be distributed in supplementary information
'''
import os
import base64

from bs4 import BeautifulSoup
import requests


ROOT = os.path.join(os.path.dirname(__file__), '..')
JSMOL = 'https://chemapps.stolaf.edu/jmol/jsmol-2019-10-30/JSmol.min.nojq.js'


def get_html():
    with open(os.path.join(ROOT, 'app', 'index.html')) as fd:
        html = fd.read()
    return BeautifulSoup(html, 'html.parser')


def get_url(url):
    if url.startswith('http://') or url.startswith('https://'):
        data = requests.get(url)
        return data.text
    else:
        if url == 'chemiscope.min.js':
            path = os.path.join(ROOT, 'dist')
        else:
            path = os.path.join(ROOT, 'app')

        with open(os.path.join(path, url)) as fd:
            return fd.read()


def inline_ressources(html):
    for script in html.find_all('script'):
        if script.has_attr('src'):
            url = script['src']
            if url.startswith('https://') or url.startswith('http://'):
                continue

            if 'JSmol.min.nojq.js' in url:
                new_script = html.new_tag('script', src=JSMOL)
                script.replace_with(new_script)
                continue

            content = get_url(url)
            src = 'data:text/javascript;base64,{}'.format(
                base64.b64encode(content.encode('utf8')).decode('utf8')
            )

            new_script = html.new_tag('script', src=src)
            if script.has_attr('defer'):
                new_script.attrs['defer'] = None
            new_script.string = content
            script.replace_with(new_script)

    for link in html.find_all('link'):
        if link.attrs.get('rel', '') == ['stylesheet']:
            url = link['href']
            if url.startswith('https://') or url.startswith('http://'):
                continue
            content = get_url(url)
            style = html.new_tag('style')
            style.string = content
            link.replace_with(style)


def main():
    html = get_html()
    inline_ressources(html)

    with open('standalone.html', 'w') as fd:
        fd.write(str(html))
        fd.write('<script id=standalone-json-data type="application/json">')

    print('created standalone.html!')


if __name__ == '__main__':
    main()
