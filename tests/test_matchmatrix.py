from unittest import TestCase
from src.matrices import build_matrix
from Bio.Seq import Seq
from numpy.testing import assert_equal


def test_should_detect_palindrome():
    sequence = Seq("ACCTAGGT")
    matrix = build_matrix(sequence, sequence, 20)
    assert_equal(matrix[:,:, 0], matrix[:,:, 1])

    
def test_should_detect_direct_repeat():
    sequence = Seq("TTACGTTACG")
    matrix = build_matrix(sequence, sequence, 20)
    assert_equal(matrix[:,:, 0], matrix[:,:, 0].T)
    assert (matrix[:,:, 0] == matrix[:,:, 0].T).all() == True
    assert (matrix[:,:, 1] == matrix[:,:, 1].T).all() == False

        
def test_should_detect_repeat():
    sequence = Seq("TTACGGCATT")
    matrix = build_matrix(sequence, sequence, 20)
    assert_equal(matrix[:,:, 0], matrix[:,:, 0].T)
    assert_equal(matrix[:,:, 1], matrix[:,:, 1].T)
    assert_equal(matrix[:,:, 2], matrix[:,:, 2].T)


class TestBuildSimpleMatrix(TestCase):
    def setUp(self):
        self.sequence = Seq("ACAT")

    def test_should_build_red_layer(self):
        layer = [
            [
                20, 0, 20, 0
            ],
            [
                0, 20, 0, 0
            ],
            [
                20, 0, 20, 0
            ],
            [
                0, 0, 0, 20
            ]
        ]
        matrix = build_matrix(self.sequence, self.sequence, 20)
        assert_equal(matrix[:,:, 0], layer)

    def test_should_build_green_layer(self):
        layer = [
            [
                20, 0, 20, 0
            ],
            [
                0, 0, 0, 20
            ], 
            [
                0, 0, 0, 0
            ],
            [
                0, 0, 0, 20
            ]
        ]

        matrix = build_matrix(self.sequence, self.sequence, 20)
        assert_equal(matrix[:,:, 1], layer)

    def test_should_build_blue_layer(self):
        layer = [

            [
                0, 20, 0, 20
            ], 
            [
                20, 0, 20, 0
            ], 
            [
                0, 20, 0, 20
            ], 
            [
                20, 20, 20, 0
            ]
        ]

        matrix = build_matrix(self.sequence, self.sequence, 20)
        assert_equal(matrix[:,:, 2], layer) 