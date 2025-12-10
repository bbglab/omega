"""Global utilities for Omega."""
import logging
import os
import subprocess
from datetime import datetime
from functools import wraps

import click
import daiquiri
from daiquiri.formatter import ColorFormatter

from omega import __logger_name__

logger = daiquiri.getLogger(__logger_name__)

DATE = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
FORMAT = '%(asctime)s - %(color)s%(levelname)-7s%(color_stop)s | %(name)s - %(color)s%(message)s%(color_stop)s'


# =========
#  Logging
# =========


def setup_logging_decorator(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        log_dir = os.path.join('.', 'log')
        command_name = click.get_current_context().command.name

        if command_name == 'run':
            cohort = click.get_current_context().params['cohort']
            fname = f'{cohort if not "None" else command_name}_{DATE}.log'
        else:
            fname = f'{command_name}_{DATE}.log'

        os.makedirs(log_dir, exist_ok=True)

        level = logging.DEBUG  if click.get_current_context().params['verbose'] else logging.INFO

        formatter = ColorFormatter(fmt=FORMAT)

        daiquiri.setup(
            level=level,
            outputs=(
                daiquiri.output.Stream(formatter=formatter),
                daiquiri.output.File(filename=os.path.join(log_dir, fname), formatter=formatter, level=logging.DEBUG),
            ),
        )

        return func(*args, **kwargs)

    return wrapper


def startup_message(version, initializing_text):
    author = 'Ferriol Calvet & Ferran Muiños @ BBGLab'
    support_email = 'bbglab@irbbarcelona.org'
    banner_width = 70

    logger.info('#' * banner_width)
    logger.info(f'{"#" + " " * (banner_width - 2) + "#"}')
    logger.info(f'{"#" + "Welcome to omega!".center(banner_width - 2) + "#"}')
    logger.info(f'{"#" + " " * (banner_width - 2) + "#"}')
    logger.info(f'{"#" + initializing_text.center(banner_width - 2) + "#"}')
    logger.info(f'{"#" + f"Version: {version}".center(banner_width - 2) + "#"}')
    logger.info(f'{"#" + f"Author: {author}".center(banner_width - 2) + "#"}')
    logger.info(f'{"#" + f"Support: {support_email}".center(banner_width - 2) + "#"}')
    logger.info(f'{"#" + " " * (banner_width - 2) + "#"}')
    logger.info('#' * banner_width)
    logger.info('')
