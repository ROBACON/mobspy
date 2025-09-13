"""Simple tests for multiline reactions."""

import pytest
from mobspy import BaseSpecies, Simulation

def test_multiline():
    # Create species
    A, B = BaseSpecies(2)
    
    # This is the multiline format that would cause issues before the fix
    
    # turn ruff format off for this block to preserve multiline layout
    # fmt: off
    (
        A >>
        B [1]
    )
    # fmt: on
    
    # Compile to test
    MySim = Simulation(A | B)
    MySim.level = -1
    result = MySim.compile()
    
    # Assert that compilation was successful (result should not be empty)
    assert result is not None
    assert len(result) > 0

def test_ruff_style_multiline_reaction():
    """Test reaction formatted the way ruff actually formats it."""
    try:
        # Create species
        A, B, C, D = BaseSpecies(4)
        
        # This is how ruff would format a long reaction:
        # Instead of: A + B >> C + D [1]
        # Ruff formats it as:
        # fmt: off
        (
            A + B >>
            C + D [1]
        )
        # fmt: on
        
        # Try to compile - this should work with our fix
        MySim = Simulation(A | B | C | D)
        MySim.level = -1
        result = MySim.compile()
        
        print("SUCCESS: Ruff-style multiline reaction works!")
        assert result is not None
        assert len(result) > 0
        
    except Exception as e:
        print(f"FAILED: {e}")
        pytest.fail(f"Ruff-style multiline reaction failed: {e}")

def test_simple_multiline_reaction():
    """Test simple reaction split after >> operator."""
    try:
        # Create species  
        A, B = BaseSpecies(2)
        
        # Simple case: A >> B [1] becomes:
        # fmt: off
        (
            A >>
            B [1]
        )
        # fmt: on
        
        MySim = Simulation(A | B)
        MySim.level = -1
        result = MySim.compile()
        
        print("SUCCESS: Simple multiline reaction works!")
        assert result is not None
        assert len(result) > 0
        
    except Exception as e:
        print(f"FAILED: {e}")
        pytest.fail(f"Simple multiline reaction failed: {e}")

def test_complex_multiline_reaction():
    """Test complex multiline reaction with stoichiometry."""
    try:
        A, B, C, D = BaseSpecies(4)
        
        # Complex reaction: 2*A + B >> 3*C + 2*D [1] becomes:
        # fmt: off
        (
            2*A + B >>
            3*C + 2*D [1]
        )
        # fmt: on
        
        MySim = Simulation(A | B | C | D)
        MySim.level = -1
        result = MySim.compile()
        
        print("SUCCESS: Complex multiline reaction works!")
        assert result is not None
        assert len(result) > 0
        
    except Exception as e:
        print(f"FAILED: {e}")
        pytest.fail(f"Complex multiline reaction failed: {e}")

def test_normal_reaction():
    """Test that normal single-line reactions still work."""
    try:
        A, B, C = BaseSpecies(3)
        A + B >> C [1]
        
        MySim = Simulation(A | B | C)
        MySim.level = -1
        result = MySim.compile()
        
        print("SUCCESS: Normal single-line reaction works!")
        assert result is not None
        assert len(result) > 0
        
    except Exception as e:
        print(f"FAILED: {e}")
        pytest.fail(f"Normal single-line reaction failed: {e}")
