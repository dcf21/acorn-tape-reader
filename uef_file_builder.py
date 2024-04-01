# -*- coding: utf-8 -*-
# uef_file_builder.py
#
# The Python script in this file produces UEF representations of Acorn tapes.
#
# Copyright (C) 2022-2024 Dominic Ford <https://dcford.org.uk/>
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
Compile and write UEF files used by Acorn emulators such as BeebEm and JSBeeb.
"""

import gzip


class UefFileBuilder:
    """
    Class to compile and write UEF files used by Acorn emulators such as BeebEm and JSBeeb.
    """

    def __init__(self):
        """
        Class to compile and write UEF files used by Acorn emulators such as BeebEm and JSBeeb.
        """

        # First compile output into a buffer
        self.output = bytearray()

        # Baud rate (used to convert time intervals into baud units)
        self.baud_rate = 1200

        # Write file header
        self.output.extend("UEF File!\0".encode("ascii"))  # ID header
        self.output.append(5)  # UEF format major version
        self.output.append(0)  # UEF format minor version

        # Write origin block
        origin_string = "acorn-tape-reader <https://github.com/dcf21/acorn-tape-reader>\0"
        self._write_chunk_header(chunk_type=0, length=len(origin_string))
        self.output.extend(origin_string.encode("ascii"))

    def write_to_file(self, filename: str) -> None:
        """
        Write UEF file to disk.

        :param:
            Filename for output UEF file.
        :return:
            None
        """

        # Write binary file
        with gzip.open(filename, "wb") as f_out:
            f_out.write(self.output)

    def _write_chunk_header(self, chunk_type: int, length: int) -> None:
        """
        Write the header at the start of a new data chunk.

        :param chunk_type:
            Numerical chunk type
        :param length:
            Length of chunk, excluding header [bytes]
        :return:
            None
        """

        # Chunk ID
        self.output.append(chunk_type & 0xFF)
        self.output.append((chunk_type >> 8) & 0xFF)

        # Chunk length
        self.output.append(length & 0xFF)
        self.output.append((length >> 8) & 0xFF)
        self.output.append((length >> 16) & 0xFF)
        self.output.append((length >> 24) & 0xFF)

    def add_silence(self, duration: float) -> None:
        """
        Add a period of silence to the output.

        :param duration:
            Duration of silence [sec]
        :return:
            None
        """

        # Convert duration into units of baud
        duration_cycles = int(duration * (self.baud_rate * 2))
        if duration_cycles < 1:
            duration_cycles = 1
        if duration_cycles > 0xFFFF:
            duration_cycles = 0xFFFF

        # Chunk type &0112
        self._write_chunk_header(chunk_type=0x112, length=2)

        # Length of silence
        self.output.append(duration_cycles & 0xFF)
        self.output.append((duration_cycles >> 8) & 0xFF)

    def add_header_tone(self, duration: float) -> None:
        """
        Add a period of header tone to the output.

        :param duration:
            Duration of header tone [sec]
        :return:
            None
        """

        # Convert duration into units of baud
        duration_cycles = int(duration * (self.baud_rate * 2))
        if duration_cycles < 1:
            duration_cycles = 1
        if duration_cycles > 0xFFFF:
            duration_cycles = 0xFFFF

        # Chunk type &0112
        self._write_chunk_header(chunk_type=0x110, length=2)

        # Length of header tone
        self.output.append(duration_cycles & 0xFF)
        self.output.append((duration_cycles >> 8) & 0xFF)

    def add_data_chunk(self, byte_list: list) -> None:
        """
        Add a data block to the output.

        :param byte_list:
            List of bytes that comprise the data block
        :return:
            None
        """

        # Chunk type &0112
        self._write_chunk_header(chunk_type=0x100, length=len(byte_list))

        # Add each byte in turn
        for item in byte_list:
            self.output.append(item['byte'])
