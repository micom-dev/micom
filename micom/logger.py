"""configures the logger for micom."""

from loguru import logger
import sys

logger.remove()
logger.add(sys.stderr, level="WARNING")
