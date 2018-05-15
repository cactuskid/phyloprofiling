import unittest
from utils import hpputils


class HpputilsTest(unittest.TestCase):

    def test_replace_characters(self):
        self.assertEqual(hpputils.replace_characters('test.replace,character(to underscore)working:well'),
                         'test_replace_character_to_underscore_working_well')


if __name__ == "__main__":
    unittest.main()