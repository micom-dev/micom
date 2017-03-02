"""configures the logger for mico."""

import logging

logger = logging.getLogger("micom")
handler = logging.StreamHandler()
formatter = logging.Formatter(
    "%(asctime)s %(name)-10s %(levelname)-8s %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)
