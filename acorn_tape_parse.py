#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# acorn_tape_parse.py
#
# This Python script extracts binary files from WAV recordings of audio
# cassette tapes recorded by 8-bit Acorn computers, including the BBC
# Micro, Acorn Electron and BBC Master computers. It can also produce
# UEF tape images for use in emulators such as BeebEm or JSBeeb.
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
This Python script extracts binary files from WAV recordings of audio cassette tapes recorded by 8-bit Acorn
computers, including the BBC Micro, Acorn Electron and BBC Master computers. It can also produce UEF tape images for
use in emulators such as BeebEm or JSBeeb.

This script was used by the author to recover all the Acorn tapes archived on the website
<https://files.dcford.org.uk/>.

By default, this script simply converts a WAV recording into UEF format, and exports all the files to a specified
output directory, together with a textual summary of the metadata associated with each file (load address, etc.).
If a more sophisticated export is required, it is simple to call the <WavAcornFileSearch> class from an external
script to perform other actions on the files found.

Usage:

* Any bit rate is supported, but >= 44.1kHz is recommended. Both mono and stereo recordings are accepted, but stereo
is recommended, and the best channel will automatically be selected -- very often one channel is (much) less noisy
than the other.

* This script currently assumes 1200 baud, as used by almost all software. It would probably be simple to add support
for the BBC Micro's 300 baud setting, but the author has no suitable test data.

* This script does not currently support Acorn Atom tapes (which use a different header format). It would probably be
simple to adapt it to do so, but the author has no suitable test data.

References:

The author used the following webpages to figure out how to decode the waveform on Acorn tapes:
https://beebwiki.mdfs.net/Acorn_cassette_format

CRC checksum calculation: https://beebwiki.mdfs.net/CRC-16
"""

import argparse
import copy
import itertools
import logging
import os
import re
import sys

from operator import itemgetter
from typing import Dict, List, Optional, Tuple

from constants import ascii
from uef_file_builder import UefFileBuilder
from wav_file_reader import WavFileReader


class WavAcornFileSearch:
    """
    Class to extract files from WAV recordings of tapes saved by 8-bit Acorn computers (e.g. the BBC Micro).
    """

    def __init__(self, input_filename: Optional[str], relative_speed: float):
        """
        Extract files from WAV recordings of tapes saved by 8-bit Acorn computers (e.g. the BBC Micro).

        :param input_filename:
            Filename of the wav file to process
        :param relative_speed:
            Assume tape runs at given fraction of normal speed (0.95 speeds up tape by 5%)
        :return:
        """

        # Input settings
        self.input_filename: Optional[str] = input_filename

        # Open wav file
        self.wav_file = WavFileReader(input_filename=self.input_filename,
                                      min_wave_amplitude_fraction=0.03)

        # Place-holder for the list of data blocks extracted from this WAV file. This pulse list can be used to
        # generate a UEF file, after called <search_wav_file>
        self.best_block_list: List[Dict] = []

        # All the phase positions we use to count wave cycles
        # These represent counting downward zero crossings, wave peaks, upward zero crossings, and wave troughs.
        self.all_phases: Tuple = (0, 90, 180, 270)

        # List of (channel, phase) configurations
        self.all_configs: List[Tuple[int, int]] = list(itertools.product(range(self.wav_file.channels),
                                                                         self.all_phases))

        # Set search settings
        self.frequency_max_variance = 0.25  # Maximum fractional variance of tone frequencies from 1200/2400 Hz
        self.relative_speed = relative_speed  # User's estimated relative speed of tape to nominal
        self.header_frequency = 2400 * relative_speed  # Hz
        self.header_period_cycles = self.wav_file.sampling_frequency / self.header_frequency  # samples per header cycle

    def search_wav_file(self):
        """
        Main entry point for searching for files from a wav recording of an Acorn computer tape.

        :return:
            List of file objects recovered
        """

        # Build a dictionary of all the file objects we recover with each configuration
        files_recovered_by_config: Dict[int, List] = {}
        blocks_recovered_by_config: Dict[int, List] = {}

        # Search for files at each phase in turn
        for config_id, (channel, phase) in enumerate(self.all_configs):
            logging.debug("Searching channel {:d} at phase {:d}".format(channel, phase))
            self.wav_file.select_channel(channel=channel)
            self.wav_file.apply_high_pass_filter(cutoff=250)  # Apply high-pass filter
            x: tuple = self.search_for_files(phase=phase)
            files_recovered_by_config[config_id] = x[0]
            blocks_recovered_by_config[config_id] = x[1]

        # Add up total number of bytes recovered with each configuration
        bytes_by_config: Dict[int, int] = {}
        for config_id in range(len(self.all_configs)):
            bytes_recovered: int = 0
            for file in files_recovered_by_config[config_id]:
                bytes_recovered += file['byte_count_without_error']
            bytes_by_config[config_id] = bytes_recovered

        # Merge the file lists we recovered with each configuration
        sorted_config_ids: List[Tuple[int, int]] = sorted(bytes_by_config.items(), key=itemgetter(1), reverse=True)

        # Store the best list of data blocks - which we may use later to generate a UEF file using <self.write_uef_file>
        best_config_id = sorted_config_ids[0][0]
        self.best_block_list = blocks_recovered_by_config[best_config_id]

        # Build merged list of all the files we recovered with each configuration
        merged_file_list = []
        timing_margin = 0.25  # maximum allowed mismatch between time position of a file seen at different phases (sec)

        # Loop over all configurations
        for config_id in [item[0] for item in sorted_config_ids]:
            # Loop over all files recovered with each configuration
            for file in files_recovered_by_config[config_id]:
                # Reject items where we didn't even successfully recover the filename
                if file['filename'] is None:
                    continue
                # Fetch the start and end time of the file on the tape
                time_start = file['start_time']
                time_end = file['final_block']['block_end_time']

                # Check if file has already been recovered at a previous config setting
                file_matches_index = None
                action = None
                for existing_file_index, existing_file in enumerate(merged_file_list):
                    # ... to match, the filename must be the same
                    if existing_file['filename'] != file['filename']:
                        continue
                    # ... to match, the end time of the new file can't be before the start of the old file
                    if time_end < existing_file['start_time'] - timing_margin:
                        continue
                    # ... to match, the start time of the new file can't be after the end of the old file
                    if time_start > existing_file['final_block']['block_end_time'] + timing_margin:
                        continue

                    # We have a match
                    file_matches_index = existing_file_index

                    # If this file failed to load successfully, and the previous instance was OK, reject the new file
                    if existing_file['is_ok'] > file['is_ok']:
                        action = None
                        break

                    # If this file loaded successfully, and the previous instance didn't, replace previous instance
                    if file['is_ok'] > existing_file['is_ok']:
                        action = "replace"
                        break

                    # If this file recovered more bytes than the previous instance, replace previous instance
                    if file['byte_count_without_error'] >= existing_file['byte_count_without_error']:
                        action = "replace append"
                        break

                    # If this file was equally as good as the previous attempt to load it, simply increment the display
                    # of which phases it was loaded at
                    action = "append"
                    break

                # We have found a new file
                if file_matches_index is None:
                    file['config_ids'] = [config_id]
                    merged_file_list.append(file)
                # We have found a better version of an existing file
                elif action == "replace":
                    file['config_ids'] = [config_id]
                    merged_file_list[file_matches_index] = file
                # We have found an equally good version of an existing file; update list of phases where we found it
                elif action == "append":
                    merged_file_list[file_matches_index]['config_ids'].append(config_id)
                # We have found a better version of an existing file, but it still didn't fully load properly
                elif action == "replace append":
                    config_ids = merged_file_list[file_matches_index]['config_ids'] + [config_id]
                    file['config_ids'] = config_ids
                    merged_file_list[file_matches_index] = file

        # Sort list of the files we found by start time, to create chronological index of the tape
        merged_file_list.sort(key=itemgetter('start_time'))

        # Return a list of all the file objects we recovered
        return merged_file_list

    def search_for_files(self, phase: int):
        """
        Search for files in a WAV audio stream, by counting wave cycles at a particular phase position in the wave
        cycle. Different tapes load better at different phase positions, due to the differing analogue audio
        chain the signal has traversed, which can introduce phase shifts. To maximise the number of files
        recovered, it is best to try all possibilities in turn.

        :return:
            List of file objects
        """

        # Check we're searching at a valid phase position
        assert phase in self.all_phases

        # Fetch list of wave cycles in the wav audio stream, measuring the start of each cycle at one of four
        # phase positions
        if phase == 0:
            wave_cycle_times = self.wav_file.fetch_zero_crossing_times(invert_wave=False)
        elif phase == 180:
            wave_cycle_times = self.wav_file.fetch_zero_crossing_times(invert_wave=True)
        elif phase == 90:
            wave_cycle_times = self.wav_file.fetch_wave_peak_times(bracket_window=self.header_period_cycles,
                                                                   invert_wave=False)
        else:
            wave_cycle_times = self.wav_file.fetch_wave_peak_times(bracket_window=self.header_period_cycles,
                                                                   invert_wave=True)

        # Make list of pulse times and lengths
        # A pulse is defined as the time interval spanned by a single wave cycle
        pulse_list = self.wav_file.fetch_pulse_list(input_events=wave_cycle_times)

        # Assign a pulse type to each pulse - does it look like 2400 Hz or 1200 Hz, or something else?
        categorised_pulse_list = self._categorise_pulse_list(pulse_list=pulse_list)

        # Extract bits from stream of pulses
        # Ones are represented by two cycles at 2400 Hz; zeros are a single 1200 Hz cycle
        bit_list = self._parse_pulse_list(pulse_list=categorised_pulse_list)

        # Extract bytes from stream of bits
        # Each byte comprises 10 bits: 0xxxxxxxx1
        byte_list = self._parse_bit_list(bit_list=bit_list)

        # Extract data blocks (header plus up to 256 bytes) from streams of bytes
        block_list = self._create_block_list(byte_list=byte_list)

        # Display a summary of all the blocks of data we found on the tape
        # self._summarise_blocks(block_list=block_list)

        # Extract a list of the files on the tape, assembling blocks together
        file_list = self._assemble_files_from_blocks(block_list=block_list)

        # Write debugging output
        self._write_debugging(pulse_list=categorised_pulse_list, bit_list=bit_list, byte_list=byte_list)

        # Return a list of the files we recovered
        return file_list, block_list

    def _categorise_pulse_list(self, pulse_list: List):
        """
        Populate the list of pulses (i.e. wave cycles) in the audio stream, determining whether they seem to be at
        2400 Hz (two cycle of which encodes a binary 1), or at 1200 Hz (which encodes a binary 0).

        :param pulse_list:
            Input list of pulses derived from <fetch_pulse_list>
        :return:
            A list of dictionaries describing the intervals.
        """

        # Minimum duration of a pulse which we categorise as a binary 0 (i.e. 1200 Hz ish)
        min_0 = 1.0 / (self.header_frequency / 2) / (1 + self.frequency_max_variance)

        # Maximum duration of a pulse which we categorise as a binary 0 (i.e. 1200 Hz ish)
        max_0 = 1.0 / (self.header_frequency / 2) * (1 + self.frequency_max_variance)

        # Minimum duration of a pulse which we categorise as a binary 1 (i.e. 2400 Hz ish)
        min_1 = 1.0 / self.header_frequency / (1 + self.frequency_max_variance)

        # Maximum duration of a pulse which we categorise as a binary 1 (i.e. 2400 Hz ish)
        max_1 = 1.0 / self.header_frequency * (1 + self.frequency_max_variance)

        # Categorise pulses as 0s, 1s, or intermediate length
        # Intermediate length pulses can arise when frequency shifts mid-cycle
        pulse_types = {
            '1': {
                'min': min_1,
                'max': max_1
            },
            'i': {
                'min': max_1,
                'max': min_0
            },
            '0': {
                'min': min_0,
                'max': max_0
            }
        }

        # Build a histogram of the number of pulses of each type
        pulse_type_histogram = {'?': 0, '0': 0, 'i': 0, '1': 0}

        # Assign a pulse type categorisation to each wave cycle
        for index, item in enumerate(pulse_list):
            pulse_length = item['length_sec']
            pulse_type = '?'
            for candidate_pulse_type, candidate_pulse_spec in pulse_types.items():
                if candidate_pulse_spec['min'] <= pulse_length <= candidate_pulse_spec['max']:
                    pulse_type = candidate_pulse_type
                    break
            item['type'] = pulse_type
            pulse_type_histogram[pulse_type] += 1

        # Log histogram of types of pulse
        logging.debug("Pulse type histogram: {}".format(repr(pulse_type_histogram)))

        # Return pulse list
        return pulse_list

    @staticmethod
    def _parse_pulse_list(pulse_list: List):
        """
        Turn a list of categorised pulses (i.e. wave cycles), into a stream of data bits.

        Zeros are encoded by a single 1200 Hz wave cycle. Ones are encoded by two 2400 Hz wave cycles. If the signal
        switches frequency mid-cycle, then intermediate-length wave cycles may arise. In this case, we make a best
        guess as to what the bits are.

        :param pulse_list:
            List of pulses to parse
        :return:
            A dictionary describing a stream of bits
        """
        # Extract bits from stream of pulses (i.e. wave cycles)
        position = 0  # Position counter in the input stream of pulses
        bit_list = []  # List of bits we've extracted
        last_bit = '1'  # The last type of bit we extracted (helps decide what to do with intermediate-length pulses)
        sequence_length = 0  # How many identical bits have we encountered in a row?
        long_sequence = False  # If we record a very long string of 1s, we've probably found a header tone

        # Cycle through the pulses we extracted from the audio stream
        while position < len(pulse_list) - 2:
            # Process pulse sequence
            pulse_time = pulse_list[position]['time']

            # Create a string describing the next three pulses
            pulse_sequence = "{}{}{}".format(pulse_list[position]['type'],
                                             pulse_list[position + 1]['type'],
                                             pulse_list[position + 2]['type'])

            # Skip over very long sequences (more than 20 wave cycles) of identical bits
            # These are probably a header tone. They never arise in real data, as all bytes are prefixed with a zero
            # and suffixed with a one.
            if sequence_length > 20 and pulse_sequence[0] == last_bit:
                long_sequence = True
                position += 1
                continue

            # Infer most likely next bit based on previous bit and the next few wave cycles
            if long_sequence and pulse_sequence[0] in (1, 'i'):
                new_bit = '1'
                advance_by = 1
            elif pulse_sequence[0] == '0':
                new_bit = '0'
                advance_by = 1
            elif pulse_sequence[0:2] in ('11', '1i'):
                new_bit = '1'
                advance_by = 2
            elif last_bit == '0' and pulse_sequence[0] == 'i':
                new_bit = '0'
                advance_by = 1
            else:
                # Invalid pulse sequence; step through pulses until we get a valid pair
                new_bit = '?'
                advance_by = 1

            # Update counter of how many identical bits we have found in succession
            if new_bit == last_bit:
                sequence_length += 1
            else:
                sequence_length = 0
                long_sequence = False
                last_bit = new_bit

            # Add new bit to the output bit list (including the time where we found it on the tape)
            bit_list.append({
                'time': pulse_time,
                'bit': new_bit
            })

            # Move to next sample
            position += advance_by

        # Output list of the bits we found
        return bit_list

    @staticmethod
    def _parse_bit_list(bit_list: List):
        """
        Take a list of the binary bits we found on the tape, into decode them into a stream of bytes. We look for
        strings of eight bits, with a zero beforehand and a one afterwards: <0xxxxxxxx1>. This is how Acorn
        computers encode each byte on the tape.

        :param bit_list:
            List of bits to assemble into bytes
        :return:
            A list of blocks, each comprising a list of bytes, separated by periods when byte-boundary synchronisation
            was lost.
        """
        # Extract bytes from stream of bits
        position = 0  # Position counter in the input stream of bits recovered from the audio stream
        synchronised = False  # Have we synchronised to a stream of bytes encoded as <0xxxxxxxx1>?
        bit_count = len(bit_list)  # How many binary bits did we find in the audio stream?
        byte_list: List[List[Dict]] = []  # A list of blocks of bytes. We start a new block whenever we lose sync.

        # Now turn stream of bits into a stream of bytes
        while position < bit_count - 9:
            # Are the next two bits of the form <0xxxxxxxx1>?
            if ((bit_list[position]['bit'] == '0') and
                    (bit_list[position + 9]['bit'] == '1') and
                    ('?' not in [i['bit'] for i in bit_list[position + 1:position + 9]])
            ):
                # If yes, then extract the value of this byte, in the range 0-255
                byte_value = sum([int(i['bit']) * pow(2, c) for c, i in enumerate(bit_list[position + 1:position + 9])])

                # If we weren't previously synchronised, then we are now! Start a new block
                if not synchronised:
                    byte_list.append([])
                # Add this byte to the current block of bytes
                byte_list[-1].append({
                    'time': bit_list[position]['time'],
                    'byte': byte_value
                })
                synchronised = True
                # Advance by 10 bits
                position += 10
            else:
                # If we didn't get a valid byte, then we've lost synchronisation
                synchronised = False
                position += 1

        # Output list of blocks of bytes
        return byte_list

    @staticmethod
    def _create_block_list(byte_list: List):
        """
        Turn a stream of bytes recovered from the audio stream into a list of Acorn data blocks. An Acorn data block
        contains a header and up to 256 bytes of data.

        :param byte_list:
            Stream of bytes to extract block information from
        :return:
            A dictionary describing the Acorn data blocks that we recovered
        """
        # Create list of Acorn data blocks
        block_list = []

        # Loop over each contiguous block of bytes we recovered
        for block_bytes in byte_list:
            # Check this block has at least some bytes
            if len(block_bytes) < 3:
                continue

            # Create template data structure to describe the Acorn data block we have found
            block_info = {
                'block_bytes': block_bytes,  # List of bytes of data, each described by a dictionary
                'block_start_time': block_bytes[0]['time'],  # Start time of the block on the tape (sec)
                'block_end_time': block_bytes[-1]['time'],  # End time of the block on the tape (sec)
                'header_pos_start': 0,  # The byte position in the block where the block header starts
                'header_pos_load_addr': 0,  # The byte position in the block where the load address was found
                'header_pos_end': 0,  # The byte position in the block where the block header ends
                'data_pos_start': 0,  # The byte position in the block where the data payload starts
                'data_pos_end': 0,  # The byte position in the block where the data payload ends
                'filename': '',  # The filename given in the block header
                'load_addr': 0,  # The 32-bit load address given in the block header
                'exec_addr': 0,  # The 32-bit execution address given in the block header
                'block_number': 0,  # The 16-bit block number given in the block header
                'data_length': 0,  # The 16-bit number of bytes in the block
                'block_flag': 0,  # The 8-bit block flags, as given in the header.
                'next_addr': 0,  # The 32-bit address of the next block, given in the block header
                'header_crc': 0,  # The header CRC (i.e. checksum) as stated on the tape
                'data_payload': [],  # The string of bytes that comprise the data payload
                'data_crc': 0,  # The data CRC (i.e. checksum) as stated on the tape
                'error': ''  # String message indicating any errors we encountered
            }
            block_list.append(block_info)

            # Check header has enough bytes, and starts with the magic byte &2A
            if block_bytes[0]['byte'] != 0x2A:
                block_info['error'] = "Wrong start byte"
                continue
            if len(block_bytes) < 21:
                block_info['error'] = "Header truncated ({:d} bytes)".format(len(block_bytes))
                continue

            # Extract filename from the header
            for char_byte in [block_bytes[i]['byte'] for i in range(1, 12)]:
                if char_byte == 0:
                    break
                block_info['filename'] += ascii[char_byte]
            if len(block_info['filename']) > 10:
                block_info['error'] = "Filename too long ({:d} chars)".format(len(block_info['filename']))
                continue

            # Check header has enough bytes
            if len(block_bytes) < 21 + len(block_info['filename']):
                block_info['error'] = "Header truncated ({:d} bytes)".format(len(block_bytes))
                continue
            block_info['header_pos_start'] = 1
            block_info['header_pos_load_addr'] = 2 + len(block_info['filename'])
            block_info['header_pos_end'] = 21 + len(block_info['filename'])

            # Read all the fields from the header
            # This assumes the header format used by the BBC Micro, also used by the Acorn Electron and BBC Master
            # If you want to read Acorn Atom tapes, you'll need to change these lines
            l = block_info['header_pos_load_addr']
            read_int = WavAcornFileSearch._read_int_from_bytes
            block_info['load_addr'] = read_int(input=block_bytes, start=l, byte_count=4)
            block_info['exec_addr'] = read_int(input=block_bytes, start=l + 4, byte_count=4)
            block_info['block_number'] = read_int(input=block_bytes, start=l + 8, byte_count=2)
            block_info['data_length'] = read_int(input=block_bytes, start=l + 10, byte_count=2)
            block_info['block_flag'] = read_int(input=block_bytes, start=l + 12, byte_count=1)
            block_info['next_addr'] = read_int(input=block_bytes, start=l + 13, byte_count=4)
            block_info['header_crc'] = read_int(input=block_bytes, start=l + 17, byte_count=2)

            # Check CRC (i.e. 16-bit checksum) of the block header
            calculated_header_crc = WavAcornFileSearch._calculate_crc(input=block_bytes,
                                                                      start=block_info['header_pos_start'],
                                                                      end=block_info['header_pos_end'] - 2)
            if calculated_header_crc != block_info['header_crc']:
                block_info['error'] = "Header CRC fail: computed {:04X}; tape says {:04X}".format(
                    calculated_header_crc, block_info['header_crc'])
                continue

            # Check block length matches the number of bytes the header told us to expect
            block_info['data_pos_start'] = block_info['header_pos_end']
            block_info['data_pos_end'] = block_info['data_pos_start'] + block_info['data_length']
            block_expected_length = block_info['data_pos_end'] + 2
            if len(block_bytes) != block_expected_length:
                block_info['error'] = "Expected {:04X} bytes; got {:04X} bytes".format(block_expected_length,
                                                                                       len(block_bytes))
                continue

            # Add data payload to the dictionary describing this Acorn data block
            block_info['data_payload'] = block_bytes[block_info['data_pos_start']:block_info['data_pos_end']]

            # Check CRC (i.e. 16-bit checksum) of the data payload
            block_info['data_crc'] = read_int(input=block_bytes,
                                              start=block_info['data_pos_end'],
                                              byte_count=2)
            calculated_data_crc = WavAcornFileSearch._calculate_crc(input=block_bytes,
                                                                    start=block_info['data_pos_start'],
                                                                    end=block_info['data_pos_end'])
            if calculated_data_crc != block_info['data_crc']:
                block_info['error'] = "Data CRC fail: computed {:04X}; tape says {:04X}".format(
                    calculated_data_crc, block_info['data_crc'])
                continue

        # Output list of Acorn data blocks
        return block_list

    @staticmethod
    def _read_int_from_bytes(input, start, byte_count):
        """
        Read an unsigned integer from a stream of bytes, with the least significant byte stored first.

        :param input:
            The input stream of bytes (each being stored as a dictionary, with the byte value stored with the key
            'byte')
        :param start:
            The position of the first byte of the integer in the byte stream
        :param byte_count:
            The number of bytes comprising the integer to be read (we can read unsigned ints of arbitrary bit widths)
        :return:
            Unsigned integer value
        """

        return sum([i['byte'] * pow(256, c) for c, i in
                    enumerate(input[start:start + byte_count])])

    @staticmethod
    def _calculate_crc(input, start, end):
        """
        Calculate a 16-bit CRC (i.e. checksum used on Acorn tapes) for a string of bytes.

        :param input:
            The input stream of bytes
        :param start:
            The position of the first byte to checksum
        :param end:
            The position of the last byte to checksum
        :return:
            16-bit checksum value
        """

        # Initialise the cyclic-redundancy check
        crc = 0
        poly = 0x1021

        # Loop over bytes
        for byte in input[start:end]:
            crc = crc ^ (byte['byte'] << 8)  # Fetch byte from memory, XOR into CRC top byte
            for i in range(8):  # Prepare to rotate 8 bits
                crc = crc << 1  # rotate
                if crc & 0x10000:  # bit 15 was set (now bit 16)...
                    crc = (crc ^ poly) & 0xFFFF  # XOR with XMODEM polynomic, and ensure CRC remains 16-bit value

        # Swap bytes around
        crc = ((crc & 0xFF) << 8) + ((crc & 0xFF00) >> 8)

        return crc  # Return updated CRC

    @staticmethod
    def summarise_blocks(block_list: List):
        """
        Display human-readable summary information about the Acorn data blocks we found on the tape.

        :param block_list:
            The list of block descriptors we found on the tape
        :return:
            None
        """

        for block in block_list:
            logging.info("[{:10.5f}] [{:50s}] [{:10s}] {:04X} {:04X} {:02X} {:08X} {:08X} {:08X}".format(
                block['block_start_time'], block['error'], block['filename'],
                block['block_number'], block['data_length'], block['block_flag'],
                block['load_addr'], block['exec_addr'], block['next_addr']
            ))

    @staticmethod
    def _assemble_files_from_blocks(block_list: List):
        """
        Extract a list of the files we found on the tape, assembling together all the contiguous streams of Acorn
        data blocks which have the same filename and consecutive block numbers.

        :param block_list:
            The list of Acorn block descriptors we found on the tape
        :return:
            None
        """

        # Start compiling a list of all the files we found on the tape
        output = []

        # Dictionary used to describe each file we find on the tape
        blank_file_descriptor = {
            'filename': None,
            'data': [],
            'start_time': None,
            'byte_count': 0,
            'byte_count_without_error': 0,
            'missed_blocks': [],
            'error_blocks': [],
            'first_block': None,
            'final_block': None,
            'message': '',
            'is_ok': False,
            'config_ids': []
        }

        # Start off with a null descriptor that we use to describe the file we're currently reading
        current_file_info = copy.deepcopy(blank_file_descriptor)

        # At the end of reading each file, check through the blocks we found and write a status message about
        # whether the file was read successfully, and in its entirety
        def populate_current_file_status_message():
            if current_file_info['filename'] is not None:
                # Bit 7 of the "block flag" byte should be set in the header of the final block of a file
                final_block_flag = current_file_info['final_block']['block_flag'] & 0x80

                # A file is deemed OK if the first block is 00, the last block is flagged, and no blocks are missing
                is_ok = ((not current_file_info['missed_blocks']) and (not current_file_info['error_blocks']) and
                         (current_file_info['first_block']['block_number'] == 0) and final_block_flag)

                # Start writing a status message for this file
                message = ""
                # Case 1: partial file with beginning and end missing
                if (current_file_info['first_block']['block_number'] > 0) and not final_block_flag:
                    message += " Blocks {:04X} - {:04X} only.".format(current_file_info['first_block']['block_number'],
                                                                      current_file_info['final_block']['block_number'])
                # Case 2: partial file with beginning missing
                elif (current_file_info['first_block']['block_number'] > 0):
                    message += " From block {:04X} only.".format(current_file_info['first_block']['block_number'])
                # Case 3: partial file with the end missing
                elif not final_block_flag:
                    message += " Truncated at block {:04X}.".format(current_file_info['final_block']['block_number'])
                # Case 4: Beginning and end of the file were present
                else:
                    message += " Last block {:04X} complete.".format(current_file_info['final_block']['block_number'])

                # Add report on any missing blocks
                if current_file_info['missed_blocks']:
                    message += " Missing blocks {}.".format(" ".join(["{:02X}".format(i)
                                                                      for i in current_file_info['missed_blocks']]))
                # Add a report on any blocks with read errors
                if current_file_info['error_blocks']:
                    message += " Errors in blocks {}.".format(" ".join(["{:02X}".format(i)
                                                                        for i in current_file_info['error_blocks']]))

                # Update file descriptor dictionary with status message and quality control flag
                current_file_info['message'] = message
                current_file_info['is_ok'] = is_ok

        # Start searching the tape for files
        # Iterate through each data block we found on the tape in turn
        for block in block_list:
            # If this block has a different filename from previous file, then finish reading old file and start new one
            if (block['filename'] != current_file_info['filename']) and not block['error']:
                # Finish reading old file
                populate_current_file_status_message()
                output.append(current_file_info)
                # Create a new file descriptor and populate it with new file's information
                current_file_info = copy.deepcopy(blank_file_descriptor)
                current_file_info['filename'] = block['filename']
                current_file_info['data'].extend([i['byte'] for i in block['data_payload']])
                current_file_info['start_time'] = block['block_start_time']
                current_file_info['byte_count'] = block['data_length']
                current_file_info['byte_count_without_error'] = block['data_length'] if not block['error'] else 0
                current_file_info['first_block'] = block
                current_file_info['final_block'] = block
            # We have another block of the file we're already reading
            elif block['filename'] == current_file_info['filename']:
                # Append the new bytes to the file
                current_file_info['data'].extend([i['byte'] for i in block['data_payload']])
                # Is this the next sequential block we were expecting? If not, report missing blocks
                current_file_info['missed_blocks'].extend(range(current_file_info['final_block']['block_number'] + 1,
                                                                block['block_number']))
                current_file_info['byte_count'] += block['data_length']
                if not block['error']:
                    current_file_info['byte_count_without_error'] += block['data_length']
                current_file_info['final_block'] = block
                # If this block failed its CRC check, report an error
                if block['error']:
                    current_file_info['error_blocks'].append(block['block_number'])
            # If the end-of-file block flag is set, then we've finished reading the file, so append it to the output
            if current_file_info['final_block'] is not None and (current_file_info['final_block']['block_flag'] & 0x80):
                populate_current_file_status_message()
                output.append(current_file_info)
                current_file_info = copy.deepcopy(blank_file_descriptor)
        # When we've finished the tape, if we have a partially read file buffered in <current_file_info>, add that to
        # our output now
        populate_current_file_status_message()
        output.append(current_file_info)

        # Return a list of all the files we found on the tape
        return output

    @staticmethod
    def extract_files(file_list: List, output_dir: str):
        """
        Write the contents of all the files we extracted from the tape into a user-supplied output directory.

        :param file_list:
            The list of files we found on the tape, as returned by <_assemble_files_from_blocks>
        :param output_dir:
            The directory in which to save the output files
        :return:
            None
        """

        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Extract each file in turn
        for index, item in enumerate(file_list):
            # Create a safe version of the filename, with illegal characters removed, and /s turned into \s
            # Also append an index to the start of the filename, since tapes commonly have multiple files with the
            # same name
            filename_safe = ("{:02d}_".format(index) +
                             re.sub('/', r'\\', item['filename'].encode('utf-8', errors='replace').decode('utf-8')))

            # Turn dots into "_dot_", since we can't save files called . or ..
            output_file = re.sub(r'\.', '_dot_', os.path.join(output_dir, filename_safe))

            # Create file
            with open(output_file, "wb") as file_handle:
                new_file_byte_array = bytearray(item['data'])
                file_handle.write(new_file_byte_array)

    def summarise_files(self, file_list: List):
        """
        Display summary information about the files we found on the tape.

        :param file_list:
            The list of files we found on the tape
        :return:
            None
        """

        # Start building output
        output = ""

        # Print column headings
        output += "[{:10s}] [{:10s}] [{:5s}] [{:10s}] {:4s} {:8s} {:8s} {:8s} {:10s} {:s}\n".format(
            "Start/sec", "End/sec", "Stat", "Filename", "Size", "LoadAddr", "ExecAddr", "NextAddr",
            "Phases", "Information"
        )

        # Print information about each file in turn
        for file_info in file_list:
            if file_info['filename'] is not None:
                # Make an indication of which configurations recovered this file
                config_indicator = ""
                for config_id in range(len(self.all_configs)):
                    if config_id in file_info['config_ids']:
                        config_indicator += "ABCDEFGH"[config_id]
                    else:
                        config_indicator += "-"

                # Write output line
                output += (
                    "[{:10.5f}] [{:10.5f}] [{:5s}] [{:10s}] {:04X} {:08X} {:08X} {:08X} [{:8s}] {:s}\n".format(
                        file_info['start_time'],
                        file_info['final_block']['block_end_time'],
                        " PASS" if file_info['is_ok'] else "*FAIL",
                        file_info['filename'], file_info['byte_count'],
                        file_info['first_block']['load_addr'],
                        file_info['first_block']['exec_addr'],
                        file_info['first_block']['next_addr'],
                        config_indicator,
                        file_info['message']
                    ))

        # Return output
        return output

    def write_uef_file(self, filename: str):
        """
        Write a .uef file representation of the data we found, enabling this tape to be loaded into emulators.

        :param filename:
            Filename of the binary UEF file we should write
        :return:
            None
        """

        if len(self.best_block_list) == 0:
            logging.info("Writing empty .uef file. Did you call <self.search_wav_file> first to parse the tape?")

        # Build UEF file
        uef_writer = UefFileBuilder()

        # Add data blocks
        current_time: float = 0  # seconds
        for block in self.best_block_list:
            # Add silence before block
            gap: float = max(0.2, block['block_start_time'] - current_time)
            if gap > 2.:
                silence_duration = min(4., gap - 2)
                uef_writer.add_silence(duration=silence_duration)

            # Add header tone before block
            header_duration = min(2., gap)
            uef_writer.add_header_tone(duration=header_duration)

            # Add contents of block
            uef_writer.add_data_chunk(byte_list=block['block_bytes'])

            # Update time
            current_time = block['block_end_time']

        # Write UEF file
        uef_writer.write_to_file(filename=filename)

    @staticmethod
    def _write_debugging(pulse_list: List, bit_list: List, byte_list: List):
        """
        Write debugging files to </tmp> describing the analysis of this wav file.

        :param pulse_list:
            List of pulses found in file
        :param bit_list:
            List of bits found in file
        :param byte_list:
            List of bytes found in file
        :return:
            None
        """
        # Output a list of all the raw pulses we found on the tape
        with open('/tmp/acorn_pulse_lengths.txt', 'wt') as f:
            f.write("# {:8s} {:10s} {}\n".format("Time/sec", "Length/ms", "Bit_type"))
            for item in pulse_list:
                f.write("{:10.5f} {:10.5f} {}\n".format(item['time'], item['length_sec'] * 1e3, item['type']))

        # Output a list of all the raw bits we found on the tape
        with open('/tmp/acorn_bits.txt', 'wt') as f:
            for item in bit_list:
                bit_time = item['time']
                bit_value = item['bit']
                f.write("[{:10.5f}] {}\n".format(bit_time, bit_value))

        # Output a list of all the raw bytes we found on the tape
        with open('/tmp/acorn_bytes.txt', 'wt') as f:
            for block in byte_list:
                for position, byte in enumerate(block):
                    byte_time = byte['time']
                    byte_value = byte['byte']

                    f.write("[{:10.5f}] {:02X} [{:s}] {:s}\n".format(
                        byte_time, byte_value, ascii[byte_value], "***" if position == 0 else "   "))


# Do it right away if we're run as a script
if __name__ == "__main__":
    # Read input parameters
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',
                        default="/tmp/box01_tape002a_acorn_eu_1992_april.wav",
                        type=str,
                        dest="input_filename",
                        help="Input WAV file to process")
    parser.add_argument('--output',
                        default="/tmp/computer_tape/",
                        type=str,
                        dest="output_directory",
                        help="Directory in which to put the extracted files")
    parser.add_argument('--uef',
                        default="/tmp/my_computer_tape.uef",
                        type=str,
                        dest="output_uef_file",
                        help="Filename for UEF file containing the contents of the tape")
    parser.add_argument('--debug',
                        action='store_true',
                        dest="debug",
                        help="Show full debugging output")
    parser.add_argument('--relative_speed',
                        default=1,
                        type=float,
                        dest="relative_speed",
                        help="Assume tape runs at given fraction of normal speed (0.95 speeds up tape by 5 percent)")
    parser.set_defaults(debug=False)
    args = parser.parse_args()

    # Set up a logging object
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO,
                        stream=sys.stdout,
                        format='[%(asctime)s] %(levelname)s:%(filename)s:%(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')
    logger = logging.getLogger(__name__)
    # logger.debug(__doc__.strip())

    # Open input audio file
    processor = WavAcornFileSearch(input_filename=args.input_filename,
                                   relative_speed=args.relative_speed)

    # Search for Acorn files
    file_list = processor.search_wav_file()

    # Print a summary of the files we found
    file_summary = processor.summarise_files(file_list=file_list)
    logging.info(file_summary)

    # Extract the Acorn files we found to output
    processor.extract_files(file_list=file_list, output_dir=args.output_directory)

    # Make UEF file
    if args.output_uef_file:
        processor.write_uef_file(filename=args.output_uef_file)
