#!/usr/bin/env python3
"""Numerical Z-Y-X gimbal-axis rank and positional-radix control.

Illustrative coordinate/axis results only: no hardware, astronomy, or new
physical law is measured. Positional radix 64 is not Base64 byte encoding.
"""
from __future__ import annotations

import json
import math

import numpy as np

DIGITS = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/'


def encode_radix64(n: int) -> str:
    if n < 0:
        return '-' + encode_radix64(-n)
    if n == 0:
        return '0'
    pieces = []
    while n:
        n, d = divmod(n, 64)
        pieces.append(DIGITS[d])
    return ''.join(reversed(pieces))


def decode_radix64(text: str) -> int:
    if text.startswith('-'):
        return -decode_radix64(text[1:])
    n = 0
    for character in text:
        n = 64 * n + DIGITS.index(character)
    return n


def body_rate_jacobian(roll: float, pitch: float) -> np.ndarray:
    """Columns map roll/pitch/yaw rates to body angular velocity for Z-Y-X."""
    return np.array([
        [1.0, 0.0, -math.sin(pitch)],
        [0.0, math.cos(roll), math.sin(roll) * math.cos(pitch)],
        [0.0, -math.sin(roll), math.cos(roll) * math.cos(pitch)],
    ])


def main() -> None:
    results = []
    for deg in (0, 45, 89, 90):
        jac = body_rate_jacobian(0.0, math.radians(deg))
        rank = int(np.linalg.matrix_rank(jac, tol=1e-12))
        results.append({
            'pitch_deg': deg,
            'rank': rank,
            'determinant': float(np.linalg.det(jac)),
            'condition_number': float(np.linalg.cond(jac)) if rank == 3 else None,
            'rank_tolerance': 1e-12,
        })
    assert [r['rank'] for r in results] == [3, 3, 3, 2]
    examples = []
    for value in (0, 10, 63, 64, 4096, 123456789012345678901234567890123):
        representation = encode_radix64(value)
        decoded = decode_radix64(representation)
        assert decoded == value
        examples.append({'decimal': str(value), 'radix_64': representation,
                         'decoded_decimal': str(decoded)})
    assert encode_radix64(64) == '10'
    print(json.dumps({
        'gimbal_convention': 'Z-Y-X, roll=0; body-rate Jacobian',
        'axis_rank_checks': results,
        'interpretation': 'At 90 degrees the conventional axes lose one independent instantaneous direction; coordinate choice does not remove physical alignment.',
        'radix_alphabet': DIGITS,
        'radix_examples': examples,
        'radix_limit': 'The same integer can be represented with different positional bases; a change of base alone is not a change of dynamics or experimental precision.',
    }, indent=2))


if __name__ == '__main__':
    main()
