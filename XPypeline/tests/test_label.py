"""Unit test for GravitySpy
"""

__author__ = 'Scott Coughlin <scott.coughlin@ligo.org>'

import os
import unittest2

class GravitSpyTests(unittest2.TestCase):
    """`TestCase` for the GravitySpy
    """
    def test_pickle(self):
        self.assertEqual(1, 1)


if __name__ == '__main__':
    unittest2.main()
