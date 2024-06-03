from docutils.parsers.rst import Directive

from .nodes import chemiscope


class Chemiscope(Directive):
    node_class = chemiscope
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {"filename": str, "mode": str, "gallery_dirs": str}

    def run(self):
        node = chemiscope()
        node["filename"] = self.options.get("filename")
        node["mode"] = self.options.get("mode")
        node["gallery_dirs"] = self.options.get("gallery_dirs")
        self.state.nested_parse(self.content, self.content_offset, node)
        return [node]
