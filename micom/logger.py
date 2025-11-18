"""configures the logger for micom."""

import logging
from rich.console import Console
from rich.logging import RichHandler

micom_console = Console()
formatter = logging.Formatter("%(message)s")
handler = RichHandler(level=logging.WARNING, markup=True, console=micom_console)
handler.setFormatter(formatter)

logging.basicConfig(
    level="WARNING", format="%(message)s", handlers=[handler]
)
