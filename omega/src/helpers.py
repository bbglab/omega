"""Helper functions for Omega CLI."""
import logging

import click

from omega import __logger_name__

LOG = logging.getLogger(__logger_name__)

# region omega init


@click.pass_context
def display_title_and_params(ctx: click.Context = None, title: str = 'Omega') -> None:
    """
    Display a styled title for the pipeline and prints click parameters if verbose is True.

    Parameters
    ----------
    title : str
        The title text to be displayed.
    verbose : bool
        Flag to indicate if click parameters should be printed.
    """
    LOG.info('=== %s ===', title)

    LOG.info('Input parameters:')
    for param, value in ctx.params.items():
        LOG.info('\t%s: %s', param, value)


# endregion omega init
