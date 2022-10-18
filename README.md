# acorn_tape_parse.py

This script in this repository extracts files from WAV recordings of cassettes saved by 8-bit Acorn computers, including the BBC Micro, Acorn Electron and BBC Master computers.

This script was used by the author to recover all the Acorn tapes archived on the website <https://files.dcford.org.uk/>.

By default, this script simply exports all the files to a specified output directory, though it is simple to call the <WavAcornFileSearch> class from an external script to perform other actions on the files found.

## Limitations

* This script only accepts 16-bit mono wav files as input.

* This script is quite sensitive to low-frequency noise. You may be able to recover more files if you use audio-editing software (e.g. Adobe Audition or Audacity) to pass the input audio through a ~100-Hz high-pass filter before calling this script.

* This script currently assumes 1200 baud, as used by almost all software. It would probably be simple to add support for the BBC Micro's 300 baud setting.

* This script does not currently support Acorn Atom tapes (which use a different header format), though it would probably be simple to adapt it to do so.

## Command-line syntax

### Usage:

./acorn_tape_parse.py [-h] [--input INPUT_FILENAME] [--output OUTPUT] [--debug] [--relative_speed RELATIVE_SPEED]

### Options:

|Switch                         |Meaning                                                                              |
|-------------------------------|-------------------------------------------------------------------------------------|
|-h, --help                     |show this help message and exit                                                      |
|--input INPUT                  |Input WAV file to process                                                            |
|--output OUTPUT                |Directory in which to put the extracted files                                        |
|--debug                        |Show full debugging output                                                           |
|--relative_speed RELATIVE_SPEED|Assume tape runs at given fraction of normal speed (0.95 speeds up tape by 5 percent)|


## License

This code is distributed under the Gnu General Public License V3. It is (C) Dominic Ford 2022.

## Author

Dominic Ford - <https://dcford.org.uk/>
