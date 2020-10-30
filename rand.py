#!/bin/env python
import json
import random
import sys

import click


def json_list(ctx, param, value):
    """Load a list from a JSON string.

    :param value: JSON string
    :return list: Loaded list
    """
    if type(value) is str:
        value = json.loads(value)
    if type(value) not in [list, type(None)]:
        ctx.fail("Invalid option: %s=%s, which must be list, str, or None." % (param.name, value))
    return value


@click.command(help="Generate a random number sequence.")
@click.option("-n", "--seqlen", type=click.IntRange(1), default=50, help="Sequence length.")
@click.option("-c", "--choices", callback=json_list, default=[1, 2, 3, 4, 5, 6], help="Choices.")
@click.option("-w", "--weights", callback=json_list, default=None, help="Weights for choices (default to uniform).")
@click.version_option("1.0.0")
def main(seqlen, choices, weights):
    choices = [str(c) for c in choices]
    print("".join(random.choices(choices, weights=weights, k=seqlen)))


if __name__ == "__main__":
    main(auto_envvar_prefix="RAND")  # pylint: disable=no-value-for-parameter,unexpected-keyword-arg
