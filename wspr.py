"""This module implements functions required to encode WSPR messages.
It is largely based on the G4JNT's excellent reference:
http://www.g4jnt.com/wspr_coding_process.pdf as well as SM0YSR's
wspr-tools: https://github.com/robertostling/wspr-tools"
The goal was to make an easy-to-understand, well-documented Python WSPR encoder
implementation."""

import string
from typing import List  # developed this on python 3.8

# fmt: off

# lookup dicts to convert characters to values and validate input
CALL_CHARS = {
    "0": 0,  "1": 1,  "2": 2,  "3": 3,  "4": 4,  "5": 5,  "6": 6,  "7": 7,  "8": 8,
    "9": 9,  "A": 10, "B": 11, "C": 12, "D": 13, "E": 14, "F": 15, "G": 16, "H": 17,
    "I": 18, "J": 19, "K": 20, "L": 21, "M": 22, "N": 23, "O": 24, "P": 25, "Q": 26,
    "R": 27, "S": 28, "T": 29, "U": 30, "V": 31, "W": 32, "X": 33, "Y": 34, "Z": 35,
    " ": 36,
}

LOC_CHARS = {
    "A": 0, "B": 1,  "C": 2,  "D": 3,  "E": 4,  "F": 5,  "G": 6,  "H": 7,  "I": 8,
    "J": 9, "K": 10, "L": 11, "M": 12, "N": 13, "O": 14, "P": 15, "Q": 16, "R": 17,
}

# polynomials used in convolution
POLY1 = 0xF2D05351
POLY2 = 0xE4613C47

# sync vector
SYNC = [
    1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
    0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0,
    1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1,
    1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0,
]

# fmt: on


class WSPRError(Exception):
    """Basic exception class for errors"""


def pack_callsign(callsign: str) -> int:
    """Packs a valid callsign in a 28 bit integer."""

    # callsigns can only be a max of 6 characters
    if len(callsign) > 6:
        raise WSPRError("Callsign larger than 6 characters")

    # callsigns can only contain digits, uppercase letters, or spaces
    for character in callsign:
        if character not in CALL_CHARS:
            raise WSPRError("Callsign contains invalid characters")

    # callsigns that don't have a number in the third spot get a space added
    # to the front
    if callsign[2] not in string.digits:
        callsign = " " + callsign

    # pad out the callsign to six digits with spaces
    while len(callsign) < 6:
        callsign += " "

    # straight from page 1 of http://www.g4jnt.com/wspr_coding_process.pdf
    packed_callsign: int = CALL_CHARS[callsign[0]]
    packed_callsign = packed_callsign * 36 + CALL_CHARS[callsign[1]]
    packed_callsign = packed_callsign * 10 + CALL_CHARS[callsign[2]]
    packed_callsign = packed_callsign * 27 + CALL_CHARS[callsign[3]] - 10
    packed_callsign = packed_callsign * 27 + CALL_CHARS[callsign[4]] - 10
    packed_callsign = packed_callsign * 27 + CALL_CHARS[callsign[5]] - 10

    return packed_callsign


def pack_locator_and_power(locator: str, power: str) -> int:
    """Packs a valid locator and power in a 22 bit integer."""

    # first two locator characters have to be A-R
    if locator[0] not in LOC_CHARS or locator[1] not in LOC_CHARS:
        raise WSPRError("First two locator characters must be between A-R")

    # second two locator characters have to be 0-9
    if locator[2] not in string.digits or locator[3] not in string.digits:
        raise WSPRError("Last two locator characters must be digits")

    # straight from page 2 of http://www.g4jnt.com/wspr_coding_process.pdf
    packed_locator: int = (
        (179 - 10 * LOC_CHARS[locator[0]] - int(locator[2])) * 180
        + 10 * LOC_CHARS[locator[1]]
        + int(locator[3])
    )

    pwr: int = int(power)
    if pwr < 0 or pwr > 60:
        raise WSPRError("Power must be int between 0-60")

    packed_locator_power: int = packed_locator * 128 + pwr + 64

    return packed_locator_power


def pack(callsign: str, locator: str, power: str) -> int:
    """Packs callsign, locator, power, and trailing zeros in an 81 bit integer"""

    packed_callsign: int = pack_callsign(callsign)
    packed_locator_power: int = pack_locator_and_power(locator, power)

    # Detailed in http://www.g4jnt.com/wspr_coding_process.pdf page 3
    # MSB -> LSB:
    # [ 28 bits of callsign ] [ 15 bits of locator ] [ 7 bits of power] [ 31 trailing zeros ]

    return ((packed_callsign << 22) + packed_locator_power) << 31


def parity(data: int) -> int:
    """Calculates the parity of an integer.
    Returns 1 for odd parity (odd number of 1s in the number) or 0 for even parity.
    Source: https://github.com/robertostling/wspr-tools/blob/master/encode.py
    This could likely be improved:
    https://www.geeksforgeeks.org/compute-parity-number-using-xor-table-look/"""

    if data < 0:
        raise WSPRError("Asked to calculate the parity of a negative number")

    parity_result: int = 0
    while data:
        parity_result ^= data & 1
        data = data >> 1
    return parity_result


def convolute(data: int) -> int:
    """Performs convolutional encoding for an 81 bit input int.
    Returns a 162 bit int. Based on the description on page 3 of
    http://www.g4jnt.com/wspr_coding_process.pdf"""

    # sanity checks
    if data < 0:
        raise WSPRError("Asked to convolute a negative number")
    if data.bit_length() > 81:
        raise WSPRError("Asked to convolute an number that is more than 81 bits")

    # initialize the output
    parity_bits: int = 0

    # initialize the registers
    reg0: int = 0
    reg1: int = 0

    # working left to right (MSB -> LSB) on our 81 bit packed int
    for shift in range(80, -1, -1):
        # pull a bit off data
        bit: int = (data >> shift) & 1

        # shift both registers one to the left
        reg0 <<= 1
        reg1 <<= 1

        # add the new bit to both registers
        reg0 += bit
        reg1 += bit

        # AND the registers, determine parity, and add the bits to the result
        parity_bits <<= 2
        parity_bits += parity(reg0 & POLY1) << 1
        parity_bits += parity(reg1 & POLY2)

    return parity_bits


def bit_reverse(byte: int) -> int:
    """Performs a bit reversal on a byte.
    For example 0110 1101 would become 1011 0110
    This could be improved with an XOR technique"""

    # sanity checks
    if byte > 255:
        raise WSPRError("Can't reverse non 8-bit numbers")

    reversed_byte: int = 0

    # swap bits 7 and 0
    reversed_byte += (byte & 0b10000000) >> 7
    reversed_byte += (byte & 0b00000001) << 7

    # swap bits 6 and 1
    reversed_byte += (byte & 0b01000000) >> 5
    reversed_byte += (byte & 0b00000010) << 5

    # swap bits 5 and 2
    reversed_byte += (byte & 0b00100000) >> 3
    reversed_byte += (byte & 0b00000100) << 3

    # swap bits 4 and 3
    reversed_byte += (byte & 0b00010000) >> 1
    reversed_byte += (byte & 0b00001000) << 1

    return reversed_byte


def interleave(data: int) -> List[int]:
    """Mixes up the bits to reduce damage from burst errors.
    Algorithm detailed on page 4 of http://www.g4jnt.com/wspr_coding_process.pdf
    Unlike the other functions this will return the result as a list of ints.
    This is useful as adding the sync vector is next."""

    source_index: int = 0
    destination: List[int] = [-1] * 162
    source: List[int] = [-1] * 162

    # create a list of bits from the 162 bit integer that was passed
    # MSB is at zero and LSB is at 161
    for i in range(0, 162):
        source[i] = (data >> (161 - i)) & 1

    for i in range(0, 256):
        dest_index: int = bit_reverse(i)
        if dest_index < 162:
            destination[dest_index] = source[source_index]
            source_index += 1
            if source_index == 162:
                return destination

    raise WSPRError("Unable to interleave all bits")


def sync(data: List[int]) -> List[int]:
    """Adds the sync vector to the data list, resulting in a symbols list (2 bits per
    symbol). Based on the algorithm on page 4 of
    http://www.g4jnt.com/wspr_coding_process.pdf"""

    symbols: List[int] = [0] * 162

    for i in range(0, 162):
        symbols[i] = SYNC[i] + 2 * data[i]

    return symbols

def encode(callsign: str, locator: str, power: str) -> List[int]:
    """Encodes a callsign, locator, and power into a list of WSPR symbols"""

    packed: int = pack(callsign, locator, power)
    #print(f"packed: {packed:#061b}")

    convoluted: int = convolute(packed)
    #print(f"convoluted: {convoluted:#0162b}")

    interleaved: List[int] = interleave(convoluted)
    #print(f"interleaved: {interleaved}")

    symbols: List[int] = sync(interleaved)
    #print(f"symbols: {symbols}")

    return symbols

#encode("ND6P", "DM04", "30")
