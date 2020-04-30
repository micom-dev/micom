"""Core objects for visualization."""

import functools
from micom.logger import logger
from os import path
from http.server import SimpleHTTPRequestHandler
from socketserver import TCPServer
from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader("micom", "data/templates"),
    autoescape=select_autoescape(["html"]),
)


class Visualization(object):
    """A visualization object.

    Attributes
    ----------
    folder : str
        The folder where the visualization was saved.
    data : dict
        The data used to create the Visualization.
    template : jinja2.Template
        The jinja template used to render the visualization.
    """

    def __init__(self, folder, data, template):
        self.folder = folder
        self.data = data
        self.template = env.get_template(template)

    def serve(self, port=8000):
        """Serve the visualization to view in a browser.

        Parameters
        ----------
        port : int
            The port on which to serve the webpage.

        Returns
        -------
        nothing
        """
        Handler = functools.partial(
            SimpleHTTPRequestHandler, directory=self.folder)
        with TCPServer(("", port), Handler) as httpd:
            logger.info("serving at port %d...", port)
            httpd.serve_forever()

    def save(self, **kwargs):
        """Render and and save the visualization."""
        out = path.join(self.folder, "index.html")
        self.template.stream(**kwargs).dump(out)
        for base, d in self.data.items():
            d.to_csv(path.join(self.folder, base + ".csv"), index=False)
