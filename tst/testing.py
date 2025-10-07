import unittest

from tst_reducible_representation import E


class TestReducibleRepresentation(unittest.TestCase):

    def setUp(self):
        self.n = 6  # Ringgröße


    def test_E(self):
        """Prüft, dass die Identitätsoperation jedes Orbital unverändert lässt"""
        for i in range(self.n):
            result = E(i)
            self.assertEqual(result, i, f"E-Test fehlgeschlagen für i={i}")


if __name__ == "__main__":
    unittest.main()
