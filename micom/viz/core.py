"""Core objects for visualization."""

import functools
from micom.logger import logger
from os import path, makedirs
from http.server import SimpleHTTPRequestHandler
from socketserver import TCPServer
from jinja2 import Environment, PackageLoader, select_autoescape

env = Environment(
    loader=PackageLoader("micom", "data/templates"),
    autoescape=select_autoescape(["html"]),
)


class VizServer(TCPServer):
    """Simple server that allows port reuse."""
    allow_reuse_address = True


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
        with VizServer(("", port), Handler) as httpd:
            logger.info("serving at port %d..." % port)
            print("Serving visualization at http://localhost:%d..." % port)
            try:
                httpd.serve_forever()
            except KeyboardInterrupt:
                print("\nStopping server.")
                httpd.shutdown()
                httpd.server_close()

    def save(self, **kwargs):
        """Render and and save the visualization."""
        makedirs(self.folder, exist_ok=True)
        out = path.join(self.folder, "index.html")
        self.template.stream(**kwargs).dump(out)
        for base, d in self.data.items():
            d.to_csv(path.join(self.folder, base + ".csv"), index=False)
