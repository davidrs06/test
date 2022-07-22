"""
Miscellaneous util functions
"""

from numpy import random
import string

def _generate_random_code(length=10, seed=None):
    """Generate random alphanumerical code, starting with a letter for proper sorting in file browsers.

    :param length: specifies length of random string, defaults to 10
    :type length: int
    :param seed: seed for random string generator, defaults to None
    :type seed: int, optional
    :return: random alphanumerical string with given length
    :rtype: string
    """

    random.seed(seed)
    return ''.join(random.choices(string.ascii_uppercase)+random.choices(string.ascii_uppercase+string.digits,k=length-1))
