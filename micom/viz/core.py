"""Core objects for visualization."""

from micom.logger import logger
from os import path
from jinja2 import Environment, PackageLoader, select_autoescape
import webbrowser

env = Environment(
    loader=PackageLoader("micom", "data/templates"),
    autoescape=select_autoescape(["html"]),
)


class Visualization(object):
    """A visualization object.

    Attributes
    ----------
    filename : str
        The filename of trhe saved visualization.
    data : dict
        The data used to create the Visualization.
    template : jinja2.Template
        The jinja template used to render the visualization.
    """

    def __init__(self, filename, data, template):
        self.filename = filename
        self.data = data
        self.template = env.get_template(template)

    def view(self):
        """Open the visualization in a browser.

        Parameters
        ----------
        None.

        Returns
        -------
        nothing
        """
        webbrowser.open("file://%s" % path.realpath(self.filename), new=2)

    def save(self, **kwargs):
        """Render and and save the visualization."""
        out = self.filename
        files = {k: d.to_csv(index=False) for k, d in self.data.items()}
        logger.info("Writing visualization to %s." % out)
        self.template.stream(files=files, **kwargs).dump(out)
