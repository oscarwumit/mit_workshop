"""
Unit and regression test for the molecool_mit package.
"""

# Import package, test suite, and other packages as needed
import molecool_mit
import pytest
import sys
import numpy as np

@pytest.fixture(scope="module")
def methane_molecule():
    symbols = np.array(['C', 'H', 'H', 'H', 'H'])

    coordinates = np.array([
        [1, 1, 1],
        [2.4, 1, 1],
        [-0.4, 1, 1],
        [1, 1, 2.4],
        [1, 1, -0.4]
    ])

    return symbols, coordinates

@pytest.mark.parametrize("p1, p2, p3, expected_value", [
    (np.array([np.sqrt(2)/2, np.sqrt(2)/2, 0]), np.array([0, 0, 0]), np.array([1, 0, 0]), 45), 
    (np.array([0, 0, -1]), np.array([0, 1, 0]), np.array([1, 0, 0]), 60 ), 
    (np.array([np.sqrt(3)/2, (1/2), 0]), np.array([0, 0, 0]), np.array([1, 0, 0]), 30), 
]
)
def test_calculate_angle_many(p1, p2, p3, expected_value):

    calculated_angle = molecool_mit.calculate_angle(p1, p2, p3, degrees=True)

    assert expected_value == pytest.approx(calculated_angle), f'{calculated_angle} {expected_value}'



def test_calculate_angle():

    r1 = np.array([0,0,-1])
    r2 = np.array([0,0,0])
    r3 = np.array([1,0,0])

    expected_value = 90

    calculated_value = molecool_mit.calculate_angle(r1, r2, r3, degrees=True)

    assert expected_value == calculated_value

def test_build_bond_list_failure(methane_molecule):

    symbols, coordinates = methane_molecule

    with pytest.raises(ValueError):
        bonds = molecool_mit.build_bond_list(coordinates, min_bond=-1)

@pytest.mark.skip
def test_build_bond_list():

    coordinates = np.array([
        [1,1,1],
        [2.4,1,1],
        [-0.4,1,1],
        [1,1,2.4],
        [1,1,-0.4],
    ])

    bonds = molecool_mit.build_bond_list(coordinates)
 
    assert len(bonds) == 4

    for bond_length in bonds.values():
        assert bond_length == 1.4

def test_calcualte_distance():
    """Test that the calculate distance function calcualtes what we expect."""

    r1 = np.array([0,0,0])
    r2 = np.array([0,1,0])

    expected_distance = 1

    observed_distance = molecool_mit.calculate_distance(r1, r2)

    assert expected_distance == observed_distance

def test_molecool_mit_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "molecool_mit" in sys.modules
