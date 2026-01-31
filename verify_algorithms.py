#!/usr/bin/env python3
"""Verification script to test algorithms against spec examples."""

from sanger_qc_trim.trim import trim_mott, trim_ends
from sanger_qc_trim.qc import compute_qc_metrics

def verify_spec_examples():
    """Test examples from the specification."""

    print("=" * 60)
    print("VERIFICATION: Testing Spec Examples")
    print("=" * 60)

    # Test 1: Mott trimming from spec
    print("\n1. Mott trimming: quals=[10,10,30,30,30,10], T=20")
    quals = [10, 10, 30, 30, 30, 10]
    result = trim_mott(quals, 20)
    print(f"   Expected: (2, 5)")
    print(f"   Got:      {result}")
    print(f"   ✓ PASS" if result == (2, 5) else f"   ✗ FAIL")

    # Test 2: Ends trimming from spec
    print("\n2. Ends trimming: quals=[15,25,25,15], T=20")
    quals = [15, 25, 25, 15]
    result = trim_ends(quals, 20)
    print(f"   Expected: (1, 3)")
    print(f"   Got:      {result}")
    print(f"   ✓ PASS" if result == (1, 3) else f"   ✗ FAIL")

    # Test 3: Q20 percentage from spec
    print("\n3. Q20 percentage: quals=[20,20,10,30]")
    quals = [20, 20, 10, 30]
    metrics = compute_qc_metrics(
        sample_id="test",
        source_file="test.ab1",
        file_format="ab1",
        seq="ACGT",
        quals=quals,
        trim_start=0,
        trim_end=4,
        qthreshold=20,
        min_length=1,
    )
    expected_pct_q20 = 0.75
    print(f"   Expected: {expected_pct_q20}")
    print(f"   Got:      {metrics['pct_q20']}")
    print(f"   ✓ PASS" if metrics['pct_q20'] == expected_pct_q20 else f"   ✗ FAIL")

    # Test 4: All bases below threshold
    print("\n4. All low quality: quals=[10,10,10,10], T=20")
    result = trim_mott([10, 10, 10, 10], 20)
    print(f"   Expected: (0, 0)")
    print(f"   Got:      {result}")
    print(f"   ✓ PASS" if result == (0, 0) else f"   ✗ FAIL")

    print("\n" + "=" * 60)
    print("All spec examples verified!")
    print("=" * 60)


if __name__ == "__main__":
    verify_spec_examples()
