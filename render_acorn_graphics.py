#!/usr/bin/env python3
# render_acorn_graphics.py
# -*- coding: utf-8 -*-
#
# The Python script in this file contains an ASCII lookup table.
#
# Copyright (C) 2022-2023 Dominic Ford <https://dcford.org.uk/>
#
# This code is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 3 of the License, or (at your option) any later
# version.
#
# You should have received a copy of the GNU General Public License along with
# this file; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA  02110-1301, USA

# ----------------------------------------------------------------------------


"""
Render a graphical representation of files containing Acorn screen memory.
"""

import argparse
import logging
import sys

from math import ceil, floor
from PIL import Image
from typing import Dict, Iterable, List, Tuple

# List of the BBC Micro's screen modes
acorn_screen_modes: Dict[int, Dict] = {
    0: {
        'width': 640,
        'bit_depth': 1
    },
    1: {
        'width': 320,
        'bit_depth': 2
    },
    2: {
        'width': 160,
        'bit_depth': 4
    },
    4: {
        'width': 320,
        'bit_depth': 1
    },
    5: {
        'width': 160,
        'bit_depth': 2
    }
}

# Color palettes to use when drawing graphics in each mode
palettes: Dict[int, List[Tuple[int]]] = {
    1: [(0, 0, 0), (255, 255, 255)],
    2: [(0, 0, 0), (255, 0, 0), (255, 255, 0), (255, 255, 255)],
    4: [(0, 0, 0), (255, 0, 0), (0, 255, 0), (255, 255, 0), (0, 0, 255), (255, 0, 255), (0, 255, 255), (255, 255, 255),
        (0, 0, 0), (255, 0, 0), (0, 255, 0), (255, 255, 0), (0, 0, 255), (255, 0, 255), (0, 255, 255), (255, 255, 255)]
}


def render_acorn_graphics(filename: str, output: str, screen_mode: int, offset: int = 0):
    """
    Render Acorn screen memory

    :param filename:
        The filename of the binary input file
    :param output:
        The filename of the graphical output
    :param screen_mode:
        The Acorn screen mode to simulate
    :param offset:
        The offset within the binary file to start rendering
    :return:
        None
    """

    assert screen_mode in acorn_screen_modes

    # Read input file
    with open(filename, "rb") as f:
        input_bytes = f.read()

    # Produce graphical output
    create_image_from_bytes(byte_list=input_bytes, output=output, screen_mode=screen_mode, offset=offset)


def create_image_from_bytes(byte_list: Iterable, output: str, screen_mode: int, offset: int = 0):
    """
    Create a listing of a input file

    :param byte_list:
        The bytes of the input file
    :param output:
        The filename of the graphical output
    :param screen_mode:
        The Acorn screen mode to simulate
    :param offset:
        The offset within the binary file to start rendering
    :return:
        None
    """

    byte_list = byte_list[offset:]

    assert screen_mode in acorn_screen_modes
    mode_info = acorn_screen_modes[screen_mode]
    palette = palettes[mode_info['bit_depth']]
    y_scaling = 2

    output_size_x = 640

    bytes_per_line = mode_info['width'] * mode_info['bit_depth']
    line_count = ceil(len(byte_list) / bytes_per_line)
    pixel_width = int(output_size_x / mode_info['width'])

    output_size_y = line_count * 8 * y_scaling

    # Create image
    map_image = Image.new('RGB', (output_size_x, output_size_y), (0, 0, 0))

    # Loop over input bytes
    for pos, byte in enumerate(byte_list):
        pos_y = floor(pos / bytes_per_line) * 8 + (pos % 8)
        pos_x_row = floor((pos % bytes_per_line) / 8)
        pixels_per_byte = int(8 / mode_info['bit_depth'])
        pos_x = int(pos_x_row * pixels_per_byte)
        accumulators = [0] * pixels_per_byte
        input_bits = int(byte)
        for bit_number in range(8):
            bit_value = input_bits % 2
            input_bits = floor(input_bits / 2)
            accumulators[bit_number % pixels_per_byte] += bit_value * pow(2, floor(bit_number / pixels_per_byte))
        for x_pixel in range(pixels_per_byte):
            for y_offset in range(y_scaling):
                for x_offset in range(pixel_width):
                    map_image.putpixel(((pos_x + x_pixel) * pixel_width + x_offset,
                                        pos_y * y_scaling + y_offset),
                                       palette[accumulators[pixels_per_byte - 1 - x_pixel]])

    # Save the final image
    map_image.save(output, "PNG")


# Do it right away if we're run as a script
if __name__ == "__main__":
    # Set up a logging object
    logging.basicConfig(level=logging.INFO,
                        stream=sys.stdout,
                        format='[%(asctime)s] %(levelname)s:%(filename)s:%(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.debug(__doc__.strip())

    # Read input parameters
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',
                        required=True,
                        type=str,
                        dest="input",
                        help="The input binary file")
    parser.add_argument('--output',
                        required=True,
                        type=str,
                        dest="output",
                        help="The output graphical file")
    parser.add_argument('--mode',
                        required=True,
                        type=int,
                        dest="mode",
                        help="The Acorn screen mode to simulate")
    parser.add_argument('--offset',
                        default=0,
                        type=int,
                        dest="offset",
                        help="The offset position within the file to start rendering")
    args = parser.parse_args()

    # Render graphics
    render_acorn_graphics(filename=args.input, output=args.output, screen_mode=args.mode, offset=args.offset)
