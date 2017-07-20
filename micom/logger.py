"""configures the logger for mico."""

import logging

logger = logging.getLogger("micom")
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "%(asctime)s | %(name)s | %(levelname)s | %(message)s")
handler.setFormatter(formatter)
handler.setLevel(logging.WARNING)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


def file_logger(filepath):
    """Log micom messages to file."""
    fh = logging.FileHandler(filepath)
    fh.setFormatter(formatter)
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)
